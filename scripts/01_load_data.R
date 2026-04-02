# scripts/01_load_data.R

clinical_patient <- read.table(
  "../data/raw/data_clinical_patient.txt",
  header = TRUE,
  sep = "\t",
  comment.char = "#",
  stringsAsFactors = FALSE,
  fill = TRUE,
  check.names = FALSE
)

clinical_sample <- read.table(
  "../data/raw/data_clinical_sample.txt",
  header = TRUE,
  sep = "\t",
  comment.char = "#",
  quote = "",
  stringsAsFactors = FALSE,
  fill = TRUE
)

expr <- read.table(
  "../data/raw/data_mrna_seq_v2_rsem.txt",
  header = TRUE,
  sep = "\t",
  comment.char = "#",
  quote = "",
  stringsAsFactors = FALSE,
  fill = TRUE,
  check.names = FALSE
)