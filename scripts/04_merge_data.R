# scripts/04_merge_data.R

# Merge da clin_clean e da sample_clean
clin_final <- merge(
  clin_clean,
  sample_clean,
  by = "PATIENT_ID"
)

# Garantir que existem amostras em comum
colnames(expr_log2) <- gsub("\\.", "-", colnames(expr_log2))

common_samples <- intersect(
  colnames(expr_log2),
  clin_final$SAMPLE_ID
)

# Filtrar matriz de expressão
expr_final <- expr_log2[, common_samples]

# Filtrar dataset clínico
clin_final <- clin_final[clin_final$SAMPLE_ID %in% common_samples, ]

# Alinhar a ordem das amostras
clin_final <- clin_final[match(common_samples, clin_final$SAMPLE_ID), ]

# Verificar alinhamento
cat("Número de amostras em comum:", length(common_samples), "\n")
cat("Alinhamento correto entre expressão e clínico:",
    all(colnames(expr_final) == clin_final$SAMPLE_ID), "\n")

# Dimensões finais
cat("Dimensão da matriz de expressão final:", dim(expr_final), "\n")
cat("Dimensão do dataset clínico final:", dim(clin_final), "\n")

saveRDS(expr_final, "../data/processed/expr_final.rds")
saveRDS(clin_final, "../data/processed/clin_final.rds")
