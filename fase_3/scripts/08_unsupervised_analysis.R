# scripts/08_unsupervised_analysis.R

library(Biobase)
library(limma)
library(pheatmap)  # substitui heatmap() base
library(Rtsne)     # t-SNE
library(uwot)      # UMAP
library(cluster)   # silhouette

# Partir do ExpressionSet

meta <- pData(eset)

# manter só amostras com idade disponível
keep <- !is.na(meta$AGE)
eset_unsup <- eset[, keep]

meta_unsup <- pData(eset_unsup)

# criar grupo etário para interpretação
mediana_idade <- median(meta_unsup$AGE, na.rm = TRUE)
meta_unsup$Grupo_Idade <- ifelse(
  meta_unsup$AGE < mediana_idade,
  "Jovem",
  "Idoso"
)
meta_unsup$Grupo_Idade <- factor(meta_unsup$Grupo_Idade, levels = c("Jovem", "Idoso"))

pData(eset_unsup) <- meta_unsup

cat("Mediana da idade (corte):", mediana_idade, "\n")
cat("Distribuição por grupo etário:\n")
print(table(pData(eset_unsup)$Grupo_Idade, useNA = "ifany"))

# Remover genes com variância zero

gene_var <- apply(exprs(eset_unsup), 1, var, na.rm = TRUE)
eset_var <- eset_unsup[gene_var > 0, ]

cat("Após remover genes com variância zero:", dim(eset_var), "\n")

# Selecionar genes mais variáveis

gene_var2 <- apply(exprs(eset_var), 1, var, na.rm = TRUE)

# para PCA/MDS/t-SNE/UMAP: top 500 (mais representativo que 100)
top_n_dimred <- min(500, nrow(eset_var))
top_genes_dimred <- names(sort(gene_var2, decreasing = TRUE))[1:top_n_dimred]

# para clustering/heatmap: compromisso entre legibilidade e informação
top_n_vis <- min(40, nrow(eset_var))
top_genes_vis <- names(sort(gene_var2, decreasing = TRUE))[1:top_n_vis]

expr_dimred <- exprs(eset_var[top_genes_dimred, ])
expr_vis <- exprs(eset_var[top_genes_vis, ])
meta_final <- pData(eset_var)

cat("Dimensão matriz PCA/MDS/t-SNE/UMAP:", dim(expr_dimred), "\n")
cat("Dimensão matriz clustering/heatmap:", dim(expr_vis), "\n")

# cores para grupo etário
grupo_col <- ifelse(meta_final$Grupo_Idade == "Jovem", "steelblue", "tomato")


# PCA das amostras

pca <- prcomp(t(expr_dimred), scale. = TRUE)
pca_summary <- summary(pca)$importance

var_pc1 <- round(100 * pca_summary[2, 1], 1)
var_pc2 <- round(100 * pca_summary[2, 2], 1)
var_pc3 <- round(100 * pca_summary[2, 3], 1)

cat("\nVariância explicada - PC1:", var_pc1, "% | PC2:", var_pc2, "% | PC3:", var_pc3, "%\n")

# Scree plot — quantas PCs são relevantes?
var_explicada <- pca_summary[2, 1:min(30, ncol(pca$x))]
barplot(
  var_explicada * 100,
  names.arg = paste0("PC", seq_along(var_explicada)),
  col = "steelblue",
  las = 2,
  main = "Scree plot — variância explicada por PC",
  xlab = "Componente Principal",
  ylab = "Variância explicada (%)"
)

# Variância acumulada
n_pcs_80 <- min(which(pca_summary[3, ] >= 0.80))
cat("Nº de PCs para explicar >= 80% da variância:", n_pcs_80, "\n")

var_acum <- pca_summary[3, 1:min(30, ncol(pca$x))]
plot(
  seq_along(var_acum), var_acum * 100,
  type = "b", pch = 19, col = "steelblue",
  xlab = "Número de PCs",
  ylab = "Variância acumulada (%)",
  main = "Variância acumulada explicada pela PCA"
)
abline(h = 80, lty = 2, col = "red")
legend("bottomright", "80% variância", lty = 2, col = "red", bty = "n")

# Score plot PC1 vs PC2
plot(
  pca$x[, 1], pca$x[, 2],
  col = grupo_col,
  pch = 16,
  xlab = paste0("PC1 (", var_pc1, "%)"),
  ylab = paste0("PC2 (", var_pc2, "%)"),
  main = "PCA das amostras"
)
legend(
  "bottomright",
  legend = c("Jovem", "Idoso"),
  col = c("steelblue", "tomato"),
  pch = 16,
  title = "Grupo Etário",
  bty = "n"
)


# MDS das amostras

# MDS clássico com cmdscale (usa distâncias euclidianas)
dist_amostras <- dist(t(expr_dimred), method = "euclidean")
mds_classico  <- cmdscale(dist_amostras, k = 2, eig = TRUE)

var_mds1 <- round(mds_classico$eig[1] / sum(abs(mds_classico$eig)) * 100, 1)
var_mds2 <- round(mds_classico$eig[2] / sum(abs(mds_classico$eig)) * 100, 1)
cat("\nVariância MDS - Dim1:", var_mds1, "% | Dim2:", var_mds2, "%\n")

plot(
  mds_classico$points[, 1], mds_classico$points[, 2],
  col = grupo_col, pch = 16,
  xlab = paste0("Coordenada 1 (", var_mds1, "%)"),
  ylab = paste0("Coordenada 2 (", var_mds2, "%)"),
  main = "MDS clássico das amostras"
)
legend(
  "bottomright",
  legend = c("Jovem", "Idoso"),
  col = c("steelblue", "tomato"),
  pch = 16,
  title = "Grupo Etário",
  bty = "n"
)

# MDS com plotMDS do limma (Leading logFC)
mds <- plotMDS(expr_dimred, plot = FALSE)

plot(
  mds$x, mds$y,
  col = grupo_col,
  pch = 16,
  xlab = "Leading logFC dim 1",
  ylab = "Leading logFC dim 2",
  main = "MDS (Leading logFC) das amostras"
)
legend(
  "bottomright",
  legend = c("Jovem", "Idoso"),
  col = c("steelblue", "tomato"),
  pch = 16,
  title = "Grupo Etário",
  bty = "n"
)


# t-SNE

cat("\n--- t-SNE ---\n")

# t-SNE não aceita linhas duplicadas — remover
expr_t    <- t(expr_dimred)
dup_rows  <- duplicated(expr_t)
expr_t_nd <- expr_t[!dup_rows, ]
meta_tsne <- meta_final[!dup_rows, ]
grupo_col_tsne <- grupo_col[!dup_rows]

cat("Amostras únicas para t-SNE:", nrow(expr_t_nd), "\n")

set.seed(42)
res_tsne <- Rtsne(
  expr_t_nd,
  dims       = 2,
  perplexity = min(30, floor((nrow(expr_t_nd) - 1) / 3)),
  verbose    = FALSE,
  max_iter   = 1000
)

plot(
  res_tsne$Y[, 1], res_tsne$Y[, 2],
  col = grupo_col_tsne, pch = 16,
  xlab = "t-SNE 1",
  ylab = "t-SNE 2",
  main = "t-SNE das amostras"
)
legend(
  "bottomright",
  legend = c("Jovem", "Idoso"),
  col = c("steelblue", "tomato"),
  pch = 16,
  title = "Grupo Etário",
  bty = "n"
)


# UMAP

cat("\n--- UMAP ---\n")

set.seed(42)
umap_res <- umap(
  expr_t_nd,
  n_neighbors  = min(15, nrow(expr_t_nd) - 1),
  min_dist     = 0.1,
  metric       = "euclidean",
  n_components = 2,
  verbose      = FALSE
)

plot(
  umap_res[, 1], umap_res[, 2],
  col = grupo_col_tsne, pch = 16,
  xlab = "UMAP 1",
  ylab = "UMAP 2",
  main = "UMAP das amostras"
)
legend(
  "bottomright",
  legend = c("Jovem", "Idoso"),
  col = c("steelblue", "tomato"),
  pch = 16,
  title = "Grupo Etário",
  bty = "n"
)


# Clustering hierárquico dos genes

dist_genes <- dist(expr_vis, method = "euclidean")
hc_genes   <- hclust(dist_genes, method = "complete")

plot(
  hc_genes,
  main = "Clustering hierárquico dos genes (complete linkage)",
  xlab = "",
  sub  = "",
  cex  = 0.8
)


# Heatmap com anotações clínicas (pheatmap)

vars_annot <- c("Grupo_Idade")
if ("STAGE_GROUP" %in% colnames(meta_final)) vars_annot <- c(vars_annot, "STAGE_GROUP")
if ("SEX"         %in% colnames(meta_final)) vars_annot <- c(vars_annot, "SEX")

annotation_col <- meta_final[, vars_annot, drop = FALSE]
rownames(annotation_col) <- colnames(expr_vis)

pheatmap(
  expr_vis,
  scale             = "row",
  show_colnames     = FALSE,
  annotation_col    = annotation_col,
  clustering_method = "complete",
  main              = "Heatmap dos genes mais variáveis (top 40)",
  fontsize_row      = 7
)


# Clustering hierárquico das amostras

dist_samples <- dist(t(expr_vis), method = "euclidean")
hc_samples   <- hclust(dist_samples, method = "complete")

plot(
  hc_samples,
  labels = FALSE,
  main   = "Clustering hierárquico das amostras",
  xlab   = "",
  sub    = ""
)

# Cortar em 2 clusters e comparar com variáveis clínicas
clusters_hier <- cutree(hc_samples, k = 2)

cat("\nClusters hierárquicos vs Grupo_Idade:\n")
print(table(clusters_hier, meta_final$Grupo_Idade))

if ("STAGE_GROUP" %in% colnames(meta_final)) {
  cat("\nClusters hierárquicos vs STAGE_GROUP:\n")
  print(table(clusters_hier, meta_final$STAGE_GROUP))
}

if ("SEX" %in% colnames(meta_final)) {
  cat("\nClusters hierárquicos vs SEX:\n")
  print(table(clusters_hier, meta_final$SEX))
}


# K-means: elbow method para escolher k

cat("\n--- K-means: elbow method ---\n")

set.seed(123)
wss <- c()
for (k in 2:10) {
  km_tmp <- kmeans(t(expr_vis), centers = k, nstart = 25)
  wss <- c(wss, sum(km_tmp$withinss))
}

plot(
  2:10, wss,
  type = "b", pch = 19, col = "steelblue",
  xlab = "Número de clusters (k)",
  ylab = "WSS total",
  main = "Elbow method — escolha do k ótimo"
)

# K-means com k=2 e nstart adequado
set.seed(123)
km <- kmeans(t(expr_vis), centers = 2, nstart = 25)

cat("\nDistribuição dos clusters:\n")
print(table(km$cluster))

cat("\nClusters vs Grupo_Idade:\n")
print(table(km$cluster, meta_final$Grupo_Idade))

if ("STAGE_GROUP" %in% colnames(meta_final)) {
  cat("\nClusters vs STAGE_GROUP:\n")
  print(table(km$cluster, meta_final$STAGE_GROUP))
}

if ("SEX" %in% colnames(meta_final)) {
  cat("\nClusters vs SEX:\n")
  print(table(km$cluster, meta_final$SEX))
}


# Guardar objetos

saveRDS(eset_unsup, "../data/processed/eset_unsup.rds")
saveRDS(expr_dimred, "../data/processed/expr_dimred_top500.rds")
saveRDS(expr_vis, "../data/processed/expr_vis_top40.rds")
saveRDS(pca, "../data/processed/pca_result.rds")
saveRDS(mds_classico, "../data/processed/mds_classico.rds")
saveRDS(res_tsne, "../data/processed/tsne_result.rds")
saveRDS(umap_res, "../data/processed/umap_result.rds")
saveRDS(km, "../data/processed/kmeans_k2.rds")
