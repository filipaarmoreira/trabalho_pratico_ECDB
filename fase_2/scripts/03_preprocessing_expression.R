# scripts/03_preprocessing_expression.R
# scripts/03_preprocessing_expression.R

# Packages Bioconductor úteis nesta fase
library(limma)
library(genefilter)


# Identificar colunas de anotação e de amostras
hugo_col <- intersect(c("Hugo_Symbol", "Hugo Symbol"), colnames(expr))
entrez_col <- intersect(c("Entrez_Gene_Id", "Entrez Gene Id"), colnames(expr))

if (length(hugo_col) == 0) stop("Coluna de símbolo génico não encontrada.")
if (length(entrez_col) == 0) stop("Coluna de Entrez não encontrada.")

hugo_col <- hugo_col[1]
entrez_col <- entrez_col[1]

cat("Coluna de símbolo génico:", hugo_col, "\n")
cat("Coluna de Entrez ID:", entrez_col, "\n")

# colunas das amostras TCGA
sample_cols <- grep("^TCGA[-.]", colnames(expr), value = TRUE)
cat("Número de colunas identificadas como amostras:", length(sample_cols), "\n")

if (length(sample_cols) == 0) stop("Nenhuma coluna de amostra TCGA encontrada.")

# Construir matriz numérica de expressão
expr_mat <- data.matrix(expr[, sample_cols, drop = FALSE])

cat("Dimensão inicial da matriz numérica:", dim(expr_mat), "\n")

# Construir IDs dos genes de forma robusta
gene_symbol <- trimws(as.character(expr[[hugo_col]]))
entrez_id   <- trimws(as.character(expr[[entrez_col]]))

# usar símbolo génico quando existir; senão cair para Entrez
gene_ids <- ifelse(
  !is.na(gene_symbol) & gene_symbol != "",
  gene_symbol,
  ifelse(!is.na(entrez_id) & entrez_id != "",
         paste0("ENTREZ_", entrez_id),
         NA)
)

rownames(expr_mat) <- gene_ids

# Remover genes sem identificador

valid_genes <- !is.na(rownames(expr_mat)) & trimws(rownames(expr_mat)) != ""
expr_mat <- expr_mat[valid_genes, , drop = FALSE]

cat("Após remover genes sem ID:", dim(expr_mat), "\n")

# Verificar e remover genes com NA

na_by_gene <- rowSums(is.na(expr_mat))
cat("Número total de NA na matriz:", sum(is.na(expr_mat)), "\n")
cat("Genes com pelo menos 1 NA:", sum(na_by_gene > 0), "\n")

# neste dataset espera-se que sejam poucos; removemos por simplicidade
expr_mat <- expr_mat[na_by_gene == 0, , drop = FALSE]

cat("Após remover genes com NA:", dim(expr_mat), "\n")

# Colapsar genes duplicados com limma

expr_mat <- limma::avereps(expr_mat, ID = rownames(expr_mat))

cat("Após colapsar genes duplicados:", dim(expr_mat), "\n")

# Filtragem de genes pouco expressos
# manter genes com expressão > 1 em pelo menos 10% das amostras
min_samples <- ceiling(0.10 * ncol(expr_mat))
keep_genes <- rowSums(expr_mat > 1) >= min_samples

expr_filt <- expr_mat[keep_genes, , drop = FALSE]

cat("Após filtrar genes pouco expressos:", dim(expr_filt), "\n")

# Transformação log para exploração

expr_log2 <- log2(expr_filt + 1)

cat("Após transformação log2:", dim(expr_log2), "\n")


# Visualizar distribuição dos SDs por gene

sds <- genefilter::rowSds(expr_log2)
m <- median(sds, na.rm = TRUE)

hist(
  sds,
  breaks = 50,
  col = "gray",
  main = "Distribuição dos desvios-padrão por gene",
  xlab = "Desvio-padrão da expressão"
)

abline(v = m, col = "red", lwd = 3, lty = 2)

legend(
  "topright",
  legend = c("Cutoff = mediana(sd)"),
  col = "red",
  lwd = 3,
  lty = 2,
  bty = "n"
)

# Filtro de flat patterns

keep_var <- sds > m
expr_log2 <- expr_log2[keep_var, , drop = FALSE]

cat("Após filtro de flat patterns (SD > mediana):", dim(expr_log2), "\n")