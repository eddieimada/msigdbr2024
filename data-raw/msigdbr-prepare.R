library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(stringr)
library(usethis)
options(pillar.print_max = 100)

# Import MSigDB gene sets -----

# Define MSigDB download variables
mdb_version <- "2024.2.Hs"
mdb_db <- str_glue("msigdb_v{mdb_version}.db")
mdb_db_zip <- str_glue("{mdb_db}.zip")
mdb_url_base <- "https://data.broadinstitute.org/gsea-msigdb/msigdb"
mdb_db_zip_url <- str_glue("{mdb_url_base}/release/{mdb_version}/{mdb_db_zip}")

# Download the MSigDB SQLite file
options(timeout = 150)
download.file(url = mdb_db_zip_url, destfile = mdb_db_zip)
unzip(mdb_db_zip, exdir = ".")
file.remove(mdb_db_zip)

# Check MSigDB SQLite file size in bytes
utils:::format.object_size(file.size(mdb_db), units = "auto")

# Open database connection to MSigDB SQLite file and extract tables as tibbles
# https://docs.gsea-msigdb.org/#MSigDB/MSigDB_SQLite_Database/
db <- DBI::dbConnect(RSQLite::SQLite(), dbname = mdb_db, flags = RSQLite::SQLITE_RO)

gene_set_source_member <- tbl(db, "gene_set_source_member") %>% as_tibble()
source_member <- tbl(db, "source_member") %>% as_tibble()
gene_set <- tbl(db, "gene_set") %>% as_tibble()
gene_symbol <- tbl(db, "gene_symbol") %>% as_tibble()
gene_set_details <- tbl(db, "gene_set_details") %>% as_tibble()
publication <- tbl(db, "publication") %>% as_tibble()

# Close database connection to MSigDB SQLite and delete the MSigDB SQLite file since it is no longer needed
DBI::dbDisconnect(db)
file.remove(mdb_db)

# Create a table for gene sets
mdb_tbl <- gene_set %>%
  inner_join(gene_set_details, join_by(id == gene_set_id)) %>%
  left_join(publication, join_by(publication_id == id)) %>%
  separate(collection_name, into = c("gs_cat", "gs_subcat"), sep = ":", remove = TRUE, extra = "merge") %>%
  select(gs_name = standard_name,
         gs_id = systematic_name,
         gs_cat,
         gs_subcat,
         gs_pmid = PMID,
         gs_geoid = GEO_id,
         gs_exact_source = exact_source,
         gs_url = external_details_URL,
         gs_description = description_brief) %>%
  filter(gs_cat != "ARCHIVED") %>%
  replace_na(list(gs_subcat = "", gs_pmid = "", gs_geoid = "", gs_exact_source = "", gs_url = "", gs_description = ""))

# Get the number of gene sets per collection (for testing)
msigdb_category_genesets <- mdb_tbl %>%
  distinct(gs_cat, gs_subcat, gs_id) %>%
  count(gs_cat, gs_subcat, name = "n_genesets")
msigdb_category_genesets

# Import MSigDB Ensembl mappings -----

# Download the MSigDB Ensembl mappings
ensembl_url <- str_glue("{mdb_url_base}/annotations/human/Human_Ensembl_Gene_ID_MSigDB.v{mdb_version}.chip")
ensembl_tbl <- read_tsv(ensembl_url, progress = FALSE, show_col_types = FALSE)
ensembl_tbl <- ensembl_tbl %>% select(human_ensembl_gene = `Probe Set ID`, human_gene_symbol = `Gene Symbol`)

# Generate a gene sets table -----

# Create a table for gene sets
msigdbr_genesets <- mdb_tbl %>%
  distinct() %>%
  arrange(gs_name, gs_id)

if (nrow(msigdbr_genesets) != sum(msigdb_category_genesets$n_genesets)) stop()

# Extract gene set members -----

# Create a table for genes in a tidy/long format (one gene per row)
geneset_genes <- gene_set %>%
  inner_join(gene_set_details, join_by(id == gene_set_id)) %>%
  select(id, systematic_name) %>%
  inner_join(gene_set_source_member, join_by(id == gene_set_id)) %>%
  inner_join(source_member, join_by(source_member_id == id)) %>%
  inner_join(gene_symbol, join_by(gene_symbol_id == id)) %>%
  select(gs_id = systematic_name, source_gene = source_id, human_entrez_gene = NCBI_id, human_gene_symbol = symbol) %>%
  replace_na(list(human_entrez_gene = 0L, human_gene_symbol = ""))
nrow(geneset_genes) %>% prettyNum(big.mark = ",")

# Check for any strange patterns
geneset_genes %>%
  count(source_gene, sort = TRUE) %>%
  head(10)
geneset_genes %>%
  count(human_gene_symbol, human_entrez_gene, sort = TRUE) %>%
  head(10)

# Get the number of members per gene set (for testing)
# Not all members map to unique genes
msigdb_geneset_members <- geneset_genes %>% count(gs_id, name = "n_members")
msigdb_geneset_members

# Confirm that gene set sizes are reasonable
if (min(msigdb_geneset_members$n_members) < 5) stop()
if (max(msigdb_geneset_members$n_members) > 3000) stop()
if (min(geneset_genes$human_entrez_gene, na.rm = TRUE) < 1) stop()

# Skip genes without an Entrez or Ensembl ID
geneset_genes <- geneset_genes %>%
  filter(human_entrez_gene > 0 | str_detect(source_gene, "^ENSG000"))
nrow(geneset_genes) %>% prettyNum(big.mark = ",")

# Keep only the relevant fields
geneset_genes <- geneset_genes %>%
  distinct(gs_id, source_gene, human_entrez_gene, human_gene_symbol)
nrow(geneset_genes) %>% prettyNum(big.mark = ",")

# Generate gene IDs -----

# Split genes based on if they include Ensembl IDs
# Starting with MSigDB 7.0, Ensembl is the platform annotation authority
# Add internal gene ID to track both Entrez and Ensembl genes
# Using Ensembl IDs as IDs for all genes resulted in a larger data file
geneset_genes_entrez <- geneset_genes %>%
  filter(str_detect(source_gene, "^ENSG000", negate = TRUE)) %>%
  distinct(gs_id, human_entrez_gene, human_gene_symbol) %>%
  mutate(gene_id = human_entrez_gene) %>%
  arrange(gs_id, gene_id)

geneset_genes_ensembl <- geneset_genes %>%
  filter(str_detect(source_gene, "^ENSG000")) %>%
  select(gs_id, human_entrez_gene, human_ensembl_gene = source_gene, human_gene_symbol) %>%
  mutate(human_gene_symbol = if_else(human_gene_symbol == "", human_ensembl_gene, human_gene_symbol)) %>%
  mutate(gene_id = str_replace(human_ensembl_gene, "ENSG000", "9")) %>%
  mutate(gene_id = as.integer(gene_id)) %>%
  arrange(gs_id, gene_id)

# Check that the gene IDs are distinct for Entrez and Ensembl tables
intersect(geneset_genes_entrez$gene_id, geneset_genes_ensembl$gene_id)

# Most gene sets should not have only some source genes as Ensembl IDs
intersect(geneset_genes_entrez$gs_id, geneset_genes_ensembl$gs_id)

# Determine unambiguous genes with only one Entrez and Ensembl ID
clean_entrez_genes <- geneset_genes_ensembl %>%
  distinct(human_entrez_gene, human_gene_symbol, human_ensembl_gene) %>%
  count(human_entrez_gene) %>%
  filter(n == 1) %>%
  pull(human_entrez_gene)
length(clean_entrez_genes)

# Use the Entrez ID for unambiguous genes
geneset_genes_ensembl <- geneset_genes_ensembl %>%
  mutate(gene_id = if_else(human_entrez_gene %in% clean_entrez_genes, human_entrez_gene, as.character(gene_id))) %>%
  arrange(gs_id, gene_id)

# Check the number of genes
nrow(geneset_genes_entrez)
n_distinct(geneset_genes_entrez$gene_id)
n_distinct(geneset_genes_entrez$human_gene_symbol)
n_distinct(geneset_genes_entrez$human_entrez_gene)
nrow(geneset_genes_ensembl)
n_distinct(geneset_genes_ensembl$gene_id)
n_distinct(geneset_genes_ensembl$human_gene_symbol)
n_distinct(geneset_genes_ensembl$human_ensembl_gene)

#### if (length(setdiff(geneset_genes_entrez$human_gene_symbol, ensembl_tbl$human_gene_symbol))) stop()

# Add Ensembl IDs to genes without them
geneset_genes_entrez <- left_join(geneset_genes_entrez, ensembl_tbl, by = "human_gene_symbol", relationship = "many-to-many")

# Check gene numbers
nrow(geneset_genes_entrez)
n_distinct(geneset_genes_entrez$human_entrez_gene)
n_distinct(geneset_genes_entrez$human_gene_symbol)
n_distinct(geneset_genes_entrez$human_ensembl_gene)

# Generate a gene set members table -----

# Combine Entrez and Ensembl genes into a single table
msigdbr_geneset_genes <-
  bind_rows(geneset_genes_entrez, geneset_genes_ensembl) %>%
  distinct(gs_id, gene_id) %>%
  arrange(gs_id, gene_id)

# Check the total number of gene set members
nrow(msigdbr_geneset_genes) %>% prettyNum(big.mark = ",")

# Check that all the original gene sets are present
if (length(setdiff(msigdb_geneset_members$gs_id, msigdbr_geneset_genes$gs_id)) > 0) stop()

# Check that most of the original gene set members converted to genes
if (nrow(msigdbr_geneset_genes) < (sum(msigdb_geneset_members$n_members) * 0.85)) stop()
genes_members_ratio = full_join(msigdb_geneset_members, count(msigdbr_geneset_genes, gs_id, name = "n_genes"), by = "gs_id")
genes_members_ratio$ratio = genes_members_ratio$n_genes / genes_members_ratio$n_members
if (min(genes_members_ratio$n_genes) < 5) stop()
if (max(genes_members_ratio$n_genes) > 2300) stop()
if (max(genes_members_ratio$ratio) > 1) stop()
if (quantile(genes_members_ratio$ratio, 0.001) < 0.3) stop()
if (quantile(genes_members_ratio$ratio, 0.1) < 0.7) stop()
if (quantile(genes_members_ratio$ratio, 0.2) < 0.9) stop()
if (quantile(genes_members_ratio$ratio, 0.3) < 0.99) stop()

# Generate a genes table -----

# Extract the unique genes
msigdbr_genes <-
  bind_rows(geneset_genes_entrez, geneset_genes_ensembl) %>%
  select(gene_id, human_gene_symbol, human_entrez_gene, human_ensembl_gene) %>%
  distinct() %>%
  arrange(human_gene_symbol, gene_id)

# Check the total number of genes
nrow(msigdbr_genes)

# Prepare package -----

# Check the size of final tables
format(object.size(msigdbr_genesets), units = "Mb")
format(object.size(msigdbr_geneset_genes), units = "Mb")
format(object.size(msigdbr_genes), units = "Mb")

# Create package data
use_data(
  msigdbr_genesets,
  msigdbr_geneset_genes,
  msigdbr_genes,
  internal = TRUE,
  overwrite = TRUE,
  compress = "xz"
)
