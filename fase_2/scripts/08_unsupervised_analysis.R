# scripts/08_unsupervised_analysis.R

library(Biobase)
library(limma)

# Partir do ExpressionSet

meta <- pData(eset)

# manter só amostras com idade disponível
keep <- !is.na(meta$AGE)
eset_unsup <- eset[, keep]

meta_unsup <- pData(eset_unsup)

# criar grupo etário para interpretação
meta_unsup$Grupo_Idade <- ifelse(
  meta_unsup$AGE < median(meta_unsup$AGE, na.rm = TRUE),
  "Jovem",
  "Idoso"
)
meta_unsup$Grupo_Idade <- factor(meta_unsup$Grupo_Idade, levels = c("Jovem", "Idoso"))

pData(eset_unsup) <- meta_unsup

cat("Distribuição por grupo etário:\n")
print(table(pData(eset_unsup)$Grupo_Idade, useNA = "ifany"))

# Remover genes com variância zero

gene_var <- apply(exprs(eset_unsup), 1, var, na.rm = TRUE)
eset_var <- eset_unsup[gene_var > 0, ]

cat("Após remover genes com variância zero:", dim(eset_var), "\n")

# Selecionar genes mais variáveis


gene_var2 <- apply(exprs(eset_var), 1, var, na.rm = TRUE)

# para PCA/MDS: mais genes
top_n_dimred <- min(100, nrow(eset_var))
top_genes_dimred <- names(sort(gene_var2, decreasing = TRUE))[1:top_n_dimred]

# para clustering/heatmap: compromisso entre legibilidade e informação
top_n_vis <- min(40, nrow(eset_var))
top_genes_vis <- names(sort(gene_var2, decreasing = TRUE))[1:top_n_vis]

expr_dimred <- exprs(eset_var[top_genes_dimred, ])
expr_vis <- exprs(eset_var[top_genes_vis, ])
meta_final <- pData(eset_var)

cat("Dimensão matriz PCA/MDS:", dim(expr_dimred), "\n")
cat("Dimensão matriz clustering/heatmap:", dim(expr_vis), "\n")

# cores para grupo etário
grupo_col <- ifelse(meta_final$Grupo_Idade == "Jovem", "black", "deeppink2")

# PCA das amostras

pca <- prcomp(t(expr_dimred), scale. = TRUE)

plot(
  pca$x[, 1], pca$x[, 2],
  col = grupo_col,
  pch = 16,
  xlab = paste0("PC1 (", round(100 * summary(pca)$importance[2, 1], 1), "%)"),
  ylab = paste0("PC2 (", round(100 * summary(pca)$importance[2, 2], 1), "%)"),
  main = "PCA das amostras"
)

legend(
  "bottomright",
  legend = c("Jovem", "Idoso"),
  col = c("black", "deeppink2"),
  pch = 16,
  title = "Grupo Etário"
)

# MDS das amostras

mds <- plotMDS(expr_dimred, plot = FALSE)

plot(
  mds$x, mds$y,
  col = grupo_col,
  pch = 16,
  xlab = "Leading logFC dim 1",
  ylab = "Leading logFC dim 2",
  main = "MDS das amostras"
)

legend(
  "bottomright",
  legend = c("Jovem", "Idoso"),
  col = c("black", "deeppink2"),
  pch = 16,
  title = "Grupo Etário"
)


# Clustering hierárquico dos genes

dist_genes <- dist(expr_vis, method = "euclidean")
hc_genes <- hclust(dist_genes, method = "complete")

plot(
  hc_genes,
  main = "Clustering hierárquico dos genes",
  xlab = "",
  sub = "",
  cex = 0.8
)

# Heatmap simples dos genes selecionados

heatmap(
  expr_vis,
  labCol = FALSE,
  cexRow = 0.8,
  main = "Heatmap dos genes mais variáveis"
)

# Clustering hierárquico das amostras

dist_samples <- dist(t(expr_vis), method = "euclidean")
hc_samples <- hclust(dist_samples, method = "complete")

plot(
  hc_samples,
  labels = FALSE,
  main = "Clustering hierárquico das amostras",
  xlab = "",
  sub = ""
)

# K-means das amostras

set.seed(123)
km <- kmeans(t(expr_vis), centers = 2)

cat("\nDistribuição dos clusters:\n")
print(table(km$cluster))

cat("\nClusters vs Grupo_Idade:\n")
print(table(km$cluster, meta_final$Grupo_Idade))


# Guardar objetos

saveRDS(eset_unsup, "../data/processed/eset_unsup.rds")
saveRDS(expr_dimred, "../data/processed/expr_dimred_top100.rds")
saveRDS(expr_vis, "../data/processed/expr_vis_top40.rds")
