# scripts/06_statistical_analysis.R

library(Biobase)
library(limma)
library(pheatmap)

# Extrair metadados do ExpressionSet

meta <- pData(eset)

# Criar grupo de estádio
meta$STAGE_GROUP <- NA

meta$STAGE_GROUP[meta$AJCC_PATHOLOGIC_TUMOR_STAGE %in% c(
  "STAGE I", "STAGE IA", "STAGE IB",
  "STAGE II", "STAGE IIA", "STAGE IIB"
)] <- "Early"

meta$STAGE_GROUP[meta$AJCC_PATHOLOGIC_TUMOR_STAGE %in% c(
  "STAGE IIIA", "STAGE IIIB", "STAGE IV"
)] <- "Advanced"

meta$STAGE_GROUP <- factor(meta$STAGE_GROUP, levels = c("Early", "Advanced"))

# Atualizar phenoData
pData(eset) <- meta

cat("Distribuição de STAGE_GROUP:\n")
print(table(pData(eset)$STAGE_GROUP, useNA = "ifany"))

# Manter só amostras com dados completos

vars_modelo <- c("STAGE_GROUP", "AGE", "SEX")
keep <- complete.cases(pData(eset)[, vars_modelo])

eset_stat <- eset[, keep]

expr_stat <- exprs(eset_stat)
meta_stat <- pData(eset_stat)

cat("\nNúmero de amostras usadas no modelo:", ncol(expr_stat), "\n")
cat("Dimensão da matriz de expressão:", dim(expr_stat), "\n")

cat("\nDistribuição final de STAGE_GROUP:\n")
print(table(meta_stat$STAGE_GROUP, useNA = "ifany"))

cat("\nDistribuição final de SEX:\n")
print(table(meta_stat$SEX, useNA = "ifany"))

summary(meta_stat$AGE)

# Modelo linear com limma

design <- model.matrix(~ STAGE_GROUP + AGE + SEX, data = meta_stat)
colnames(design) <- make.names(colnames(design))

cat("\nMatriz de desenho:\n")
print(head(design))

fit <- lmFit(expr_stat, design)
fit <- eBayes(fit)

# contraste principal: Advanced vs Early
resultados <- topTable(
  fit,
  coef = "STAGE_GROUPAdvanced",
  number = Inf,
  sort.by = "P"
)

resultados$Gene <- rownames(resultados)
resultados <- resultados[, c("Gene", setdiff(colnames(resultados), "Gene"))]

cat("\nNúmero total de genes analisados:", nrow(resultados), "\n")
cat("Genes com adj.P.Val < 0.05:",
    sum(resultados$adj.P.Val < 0.05, na.rm = TRUE), "\n")

cat("Genes com adj.P.Val < 0.01:",
    sum(resultados$adj.P.Val < 0.01, na.rm = TRUE), "\n")

# Guardar resultados

write.csv(
  resultados,
  "../data/processed/resultados_limma_stage_age_sex.csv",
  row.names = FALSE
)

write.csv(
  resultados[resultados$adj.P.Val < 0.05, ],
  "../data/processed/genes_significativos_stage_age_sex.csv",
  row.names = FALSE
)

# Volcano plot

with(resultados, plot(
  logFC,
  -log10(P.Value),
  pch = 16,
  cex = 0.6,
  main = "Volcano plot: Advanced vs Early",
  xlab = "logFC",
  ylab = "-log10(p-value)"
))

sig <- resultados$adj.P.Val < 0.05 & !is.na(resultados$adj.P.Val)

with(resultados[sig, ], points(
  logFC,
  -log10(P.Value),
  pch = 16,
  col = "red"
))

abline(h = -log10(0.05), lty = 2, col = "gray40")

# Tabela dos genes mais significativos

cat("\nTop 20 genes mais significativos:\n")
print(head(resultados, 20))

# Heatmap dos genes mais significativos

genes_sig <- resultados$Gene[resultados$adj.P.Val < 0.05]

if (length(genes_sig) >= 2) {
  
  n_heat <- min(30, length(genes_sig))
  genes_heat <- genes_sig[1:n_heat]
  
  annotation_col <- data.frame(
    STAGE_GROUP = meta_stat$STAGE_GROUP,
    AGE = meta_stat$AGE,
    SEX = meta_stat$SEX
  )
  rownames(annotation_col) <- colnames(expr_stat)
  
  pheatmap(
    expr_stat[genes_heat, , drop = FALSE],
    scale = "row",
    show_colnames = FALSE,
    annotation_col = annotation_col,
    main = "Top genes DE: Advanced vs Early"
  )
  
} else {
  cat("\nMenos de 2 genes significativos; heatmap não gerado.\n")
}

# MA-plot adaptado com limma

plotMD(
  fit,
  column = which(colnames(design) == "STAGE_GROUPAdvanced"),
  status = decideTests(fit)[, "STAGE_GROUPAdvanced"],
  main = "MD plot: Advanced vs Early"
)

