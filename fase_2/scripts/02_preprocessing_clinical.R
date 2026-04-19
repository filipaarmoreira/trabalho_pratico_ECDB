# scripts/02_preprocessing_clinical.R

# Selecionar colunas relevantes
clin_clean <- clinical_patient[, c(
  "PATIENT_ID",
  "CANCER_TYPE_ACRONYM",
  "AGE",
  "SEX",
  "AJCC_PATHOLOGIC_TUMOR_STAGE",
  "PATH_T_STAGE",
  "PATH_N_STAGE",
  "PATH_M_STAGE",
  "OS_STATUS",
  "OS_MONTHS",
  "PFS_STATUS",
  "PFS_MONTHS",
  "GENETIC_ANCESTRY_LABEL"
)]

sample_clean <- clinical_sample[, c(
  "PATIENT_ID",
  "SAMPLE_ID",
  "ANEUPLOIDY_SCORE",
  "SAMPLE_TYPE",
  "MSI_SCORE_MANTIS",
  "TMB_NONSYNONYMOUS"
)]

# Substituir strings vazias por NA
clean_empty_strings <- function(df) {
  df[] <- lapply(df, function(x) {
    if (is.character(x)) {
      x <- trimws(x)
      x[x == ""] <- NA
    }
    x
  })
  df
}

clin_clean <- clean_empty_strings(clin_clean)
sample_clean <- clean_empty_strings(sample_clean)

# Remover linhas sem identificador
clin_clean <- clin_clean[!is.na(clin_clean$PATIENT_ID), ]
sample_clean <- sample_clean[!is.na(sample_clean$PATIENT_ID) & !is.na(sample_clean$SAMPLE_ID), ]

# Converter variáveis numéricas
clin_clean$AGE <- as.numeric(clin_clean$AGE)
clin_clean$OS_MONTHS <- as.numeric(clin_clean$OS_MONTHS)
clin_clean$PFS_MONTHS <- as.numeric(clin_clean$PFS_MONTHS)
sample_clean$ANEUPLOIDY_SCORE <- as.numeric(sample_clean$ANEUPLOIDY_SCORE)
sample_clean$MSI_SCORE_MANTIS <- as.numeric(sample_clean$MSI_SCORE_MANTIS)
sample_clean$TMB_NONSYNONYMOUS <- as.numeric(sample_clean$TMB_NONSYNONYMOUS)

# Limpeza das variáveis status
clin_clean$OS_STATUS <- sub("^[0-9]+:", "", clin_clean$OS_STATUS)
clin_clean$PFS_STATUS <- sub("^[0-9]+:", "", clin_clean$PFS_STATUS)


# Converter variáveis categóricas para factor
clin_clean$CANCER_TYPE_ACRONYM <- as.factor(clin_clean$CANCER_TYPE_ACRONYM)
clin_clean$SEX <- as.factor(clin_clean$SEX)
clin_clean$AJCC_PATHOLOGIC_TUMOR_STAGE <- as.factor(clin_clean$AJCC_PATHOLOGIC_TUMOR_STAGE)
clin_clean$PATH_T_STAGE <- as.factor(clin_clean$PATH_T_STAGE)
clin_clean$PATH_N_STAGE <- as.factor(clin_clean$PATH_N_STAGE)
clin_clean$PATH_M_STAGE <- as.factor(clin_clean$PATH_M_STAGE)
clin_clean$OS_STATUS <- as.factor(clin_clean$OS_STATUS)
clin_clean$PFS_STATUS <- as.factor(clin_clean$PFS_STATUS)
clin_clean$GENETIC_ANCESTRY_LABEL <- as.factor(clin_clean$GENETIC_ANCESTRY_LABEL)
sample_clean$SAMPLE_TYPE <- as.factor(sample_clean$SAMPLE_TYPE)

# Verificações rápidas
cat("Dimensão do dataset clinical_patient limpo:", dim(clin_clean), "\n")
cat("Dimensão do dataset clinical_sample limpo:", dim(sample_clean), "\n")
summary(clin_clean$AGE)
summary(clin_clean$OS_MONTHS)
table(clin_clean$SEX, useNA = "ifany")
table(clin_clean$AJCC_PATHOLOGIC_TUMOR_STAGE, useNA = "ifany")
table(clin_clean$OS_STATUS, useNA = "ifany")
table(clin_clean$PFS_STATUS, useNA = "ifany")
table(sample_clean$SAMPLE_TYPE, useNA = "ifany")

colSums(is.na(clin_clean))
colSums(is.na(sample_clean))

sum(duplicated(clin_clean$PATIENT_ID))
sum(duplicated(sample_clean$SAMPLE_ID))

str(clin_clean)
str(sample_clean)
