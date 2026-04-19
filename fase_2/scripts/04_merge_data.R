# scripts/04_merge_data.R

library(Biobase)

# Juntar metadados clínicos
clin_final <- merge(
  sample_clean,
  clin_clean,
  by = "PATIENT_ID",
  all = FALSE,
  sort = FALSE
)

# garantir unicidade
stopifnot(!anyDuplicated(clin_final$PATIENT_ID))
stopifnot(!anyDuplicated(clin_final$SAMPLE_ID))

# Harmonizar nomes das amostras na expressão

colnames(expr_log2) <- gsub("\\.", "-", colnames(expr_log2))

# Encontrar amostras em comum

common_samples <- intersect(colnames(expr_log2), clin_final$SAMPLE_ID)

cat("Número de amostras em comum:", length(common_samples), "\n")

# Filtrar expressão e metadados

expr_final <- expr_log2[, common_samples, drop = FALSE]

clin_final <- clin_final[clin_final$SAMPLE_ID %in% common_samples, , drop = FALSE]
clin_final <- clin_final[match(common_samples, clin_final$SAMPLE_ID), , drop = FALSE]


# Verificar alinhamento

stopifnot(identical(colnames(expr_final), clin_final$SAMPLE_ID))

cat("Dimensão expr_final:", dim(expr_final), "\n")
cat("Dimensão clin_final:", dim(clin_final), "\n")
cat("Alinhamento correto entre expressão e clínico:",
    identical(colnames(expr_final), clin_final$SAMPLE_ID), "\n")

# Criar ExpressionSet

rownames(clin_final) <- clin_final$SAMPLE_ID

eset <- ExpressionSet(
  assayData = expr_final,
  phenoData = AnnotatedDataFrame(clin_final)
)

eset
