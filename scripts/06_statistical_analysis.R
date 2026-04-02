# scripts/06_statistical_analysis.R

library(limma)

# Preparar dados
amostras_com_idade <- !is.na(clin_final$AGE)
clin_stat <- clin_final[amostras_com_idade, ]
expr_stat <- expr_final[, amostras_com_idade, drop = FALSE]

cat("Número de amostras com idade disponível:", nrow(clin_stat), "\n")
cat("Dimensão da matriz de expressão para análise:", dim(expr_stat), "\n")

# Criar grupo etário
clin_stat$Grupo_Idade <- ifelse(
  clin_stat$AGE < median(clin_stat$AGE, na.rm = TRUE),
  "Jovem",
  "Idoso"
)

clin_stat$Grupo_Idade <- factor(clin_stat$Grupo_Idade, levels = c("Jovem", "Idoso"))

print(table(clin_stat$Grupo_Idade, useNA = "ifany"))


# Análise estatística univariada (exemplo com 1 gene)
gene_teste <- expr_stat[1, ]
nome_do_gene <- rownames(expr_stat)[1]

boxplot(
  gene_teste ~ clin_stat$Grupo_Idade,
  main = paste("Expressão do gene", nome_do_gene, "por grupo etário"),
  ylab = "Nível de expressão",
  xlab = "Grupo etário",
  col = c("lightgreen", "lightcoral")
)

teste_gene <- t.test(gene_teste ~ clin_stat$Grupo_Idade)
print(teste_gene)

# Expressão diferencial com limma
# (mais apropriado para dados contínuos / RSEM)
design <- model.matrix(~ Grupo_Idade, data = clin_stat)
colnames(design) <- c("Intercepto", "Idoso_vs_Jovem")

fit <- lmFit(expr_stat, design)
fit <- eBayes(fit)

resultados <- topTable(
  fit,
  coef = "Idoso_vs_Jovem",
  number = Inf,
  sort.by = "P"
)

# Adicionar coluna gene
resultados$Gene <- rownames(resultados)

# Reordenar
resultados <- resultados[, c("Gene", setdiff(colnames(resultados), "Gene"))]

cat("Número de genes analisados:", nrow(resultados), "\n")
cat("Genes com adj.P.Val < 0.05:", sum(resultados$adj.P.Val < 0.05, na.rm = TRUE), "\n")

head(resultados)

# Volcano plot simples
with(resultados, plot(
  logFC,
  -log10(P.Value),
  pch = 16,
  cex = 0.5,
  main = "Volcano plot: Idoso vs Jovem",
  xlab = "logFC",
  ylab = "-log10(p-value)"
))

sig <- resultados$adj.P.Val < 0.05 & !is.na(resultados$adj.P.Val)
with(resultados[sig, ], points(logFC, -log10(P.Value), pch = 16, col = "red"))

# Heatmap
library(pheatmap)
genes_top <- resultados$Gene[resultados$adj.P.Val < 0.05]
pheatmap(expr_log2[genes_top[1:30], ],
         scale = "row",
         show_colnames = FALSE)

# Filtrar genes significativos
genes_sig <- resultados[resultados$adj.P.Val < 0.05, ]

cat("Número de genes significativos:", nrow(genes_sig), "\n")

# Guardar lista
write.csv(genes_sig, "../data/processed/genes_significativos.csv")

# Guardar só nomes dos genes (clean)
write.table(
  rownames(genes_sig),
  "../data/processed/genes_significativos_lista.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)