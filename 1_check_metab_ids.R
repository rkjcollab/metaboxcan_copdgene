# SDS 20241217

# Check metabolite ID matching between METSIM dbs and COPDGene metabolomics.

# Setup ------------------------------------------------------------------------

setwd(paste0(Sys.getenv("RKJCOLLAB"), "/Collabs/metaboxcan"))

library(tidyverse)
library(conflicted)
library(RSQLite)

# Load data --------------------------------------------------------------------

# METSIM db
db <- dbConnect(
  RSQLite::SQLite(),
  "raw/dbs/metsim-invnorm-softimpute.db")
dbListTables(db)
# "extra"   "weights"
db_ext <- dbReadTable(db, "extra")
db_wt <- dbReadTable(db, "weights")
dbDisconnect(db)

# Additional file provided with db
db_2 <- read_delim(gzfile("raw/dbs/metsim-invnorm-softimpute.txt.gz"))
n_distinct(db_2$GENE)  # 1192
# Not sure what this file is

# COPDgene metabolomics
metab <- read_csv(
  "raw/metab/COPDGene_P2_LT20miss_knnImp_NoOut_metabs_20211021.csv")
metab_info <- read_csv(
  "raw/metab/COPDGene_P2_MetaboliteInformation_20211021.csv")


# Check IDs --------------------------------------------------------------------

# Looks like column names are metabolite IDs
metab_ids <- grep("sid", colnames(metab), invert = T, value = T)
  # these match metab_id in metabolie information file, which is based on the
  # CHEM_ID from Metabolon. CHEM_ID column also in information file

n_distinct(db_wt$gene)  # 1192 distinct gene names
n_distinct(db_ext$gene)  # 1191 gene names (not sure why one less)
setdiff(db_wt$gene, db_ext$gene)  # "C100020005"

db_ids <- unique(db_wt$gene)

length(intersect(metab_ids, db_ids))  # 0 overlap as is

# See if mapping can be done like this:
  # metabolite data: metab_id --> CHEM_ID
  # db: remove C from IDs, match to CHEM_ID
metab_info_filt <- metab_info %>%
  dplyr::filter(metab_id %in% metab_ids) %>%
  dplyr::select(metab_id, CHEM_ID)

db_ids_mod <- gsub("^C", "", db_ids)

length(intersect(metab_info_filt$CHEM_ID, db_ids_mod))  # 789
  # Based on Erika's thesis, this overlap seems correct

# How many in one but not the other?
length(setdiff(metab_info_filt$CHEM_ID, db_ids_mod))  # 206 in COPDGene, not dbs
length(setdiff(db_ids_mod, metab_info_filt$CHEM_ID))  # 403 in dbs, not in COPDGene

# See format of IDs  that don't match
metab_info_filt_no_match <- metab_info_filt %>%
  dplyr::filter(!CHEM_ID %in% db_ids_mod)
db_ids_mod_no_match <- as.data.frame(db_ids_mod) %>%
  dplyr::filter(!db_ids_mod %in% metab_info_filt$CHEM_ID)

# Doesn't appear to be because the format between the metab_id and CHEM_ID
# changes, likely just because these aren't present in the dbs
