# scripts/03_preprocessing_expression.R

# Identificar coluna com genes
gene_col <- colnames(expr)[1]

cat("Coluna de genes assumida:", gene_col, "\n")

# Guardar nomes dos genes
gene_ids <- expr[[gene_col]]

# Remover colunas de anotação extra, se existirem
sample_cols <- grep("^TCGA-", colnames(expr), value = TRUE)

cat("Número de colunas identificadas como amostras:", length(sample_cols), "\n")

expr_mat <- expr[, sample_cols, drop = FALSE]

# Converter para matriz numérica
expr_mat <- as.data.frame(lapply(expr_mat, as.numeric))
expr_mat <- as.matrix(expr_mat)

# Associar rownames
rownames(expr_mat) <- gene_ids

cat("Dimensão da matriz numérica:", dim(expr_mat), "\n")

# Remover genes sem identificador
valid_genes <- !is.na(rownames(expr_mat)) & rownames(expr_mat) != ""
expr_mat <- expr_mat[valid_genes, ]

cat("Após remover genes sem ID:", dim(expr_mat), "\n")

# Remover genes duplicados (mantém a primeira ocorrência)
expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), ]

cat("Após remover genes duplicados:", dim(expr_mat), "\n")

# Verificar missing values
cat("Número total de NA na matriz:", sum(is.na(expr_mat)), "\n")

# Remover genes com qualquer NA
expr_mat <- expr_mat[complete.cases(expr_mat), ]

cat("Após remover genes com NA:", dim(expr_mat), "\n")

# Filtragem de genes pouco expressos
min_samples <- ceiling(0.10 * ncol(expr_mat))

keep_genes <- rowSums(expr_mat > 1) >= min_samples
expr_filt <- expr_mat[keep_genes, ]

cat("Após filtrar genes pouco expressos:", dim(expr_filt), "\n")

# Transformação log para exploração
expr_log2 <- log2(expr_filt + 1)

# Verificações rápidas
cat("\nResumo dos valores da matriz filtrada:\n")
print(summary(as.vector(expr_filt)))

cat("\nResumo dos valores da matriz log2:\n")
print(summary(as.vector(expr_log2)))

# Distribuição por amostra
sample_medians <- apply(expr_log2, 2, median, na.rm = TRUE)
sample_means   <- apply(expr_log2, 2, mean, na.rm = TRUE)

cat("\nResumo das medianas por amostra:\n")
print(summary(sample_medians))

cat("\nResumo das médias por amostra:\n")
print(summary(sample_means))
