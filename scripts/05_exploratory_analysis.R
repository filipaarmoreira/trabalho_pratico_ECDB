# scripts/05_exploratory_analysis.R

# 1. Checks básicos
cat("Dimensão do dataset clínico final:", dim(clin_final), "\n")
cat("Dimensão da matriz de expressão final:", dim(expr_final), "\n")

# 2. Variáveis categóricas
tabela_estadios <- table(clin_final$AJCC_PATHOLOGIC_TUMOR_STAGE, useNA = "ifany")

barplot(
  tabela_estadios,
  main = "Distribuição dos estádios do cancro",
  col = "coral",
  las = 2,
  cex.names = 0.7
)

tabela_sexo <- table(clin_final$SEX, useNA = "ifany")
print(tabela_sexo)

barplot(
  tabela_sexo,
  main = "Distribuição por sexo",
  col = "lightpink"
)

# 2. Variáveis numéricas clínicas
summary(clin_final$AGE)

hist(
  clin_final$AGE,
  main = "Distribuição da idade dos pacientes",
  xlab = "Idade (anos)",
  ylab = "Frequência",
  col = "lightblue"
)

boxplot(
  clin_final$AGE,
  main = "Boxplot da idade",
  ylab = "Idade",
  col = "lightblue"
)

summary(clin_final$OS_MONTHS)

hist(
  clin_final$OS_MONTHS,
  main = "Distribuição de OS_MONTHS",
  xlab = "Overall survival (meses)",
  col = "lightgreen"
)

boxplot(
  clin_final$OS_MONTHS,
  main = "Boxplot de OS_MONTHS",
  ylab = "Meses",
  col = "lightgreen"
)

# 3. Variáveis moleculares
summary(clin_final$ANEUPLOIDY_SCORE)
summary(clin_final$TMB_NONSYNONYMOUS)
summary(clin_final$MSI_SCORE_MANTIS)

hist(
  clin_final$ANEUPLOIDY_SCORE,
  main = "Distribuição de Aneuploidy Score",
  xlab = "Score",
  col = "lightyellow"
)

boxplot(
  clin_final$ANEUPLOIDY_SCORE,
  main = "Boxplot - Aneuploidy Score",
  ylab = "Score",
  col = "lightyellow"
)

hist(
  clin_final$TMB_NONSYNONYMOUS,
  main = "Distribuição de TMB",
  xlab = "TMB",
  col = "lavender"
)

boxplot(
  clin_final$TMB_NONSYNONYMOUS,
  main = "Boxplot de TMB",
  ylab = "TMB",
  col = "lavender"
)

hist(
  clin_final$MSI_SCORE_MANTIS,
  main = "Distribuição de MSI_SCORE_MANTIS",
  xlab = "MSI_SCORE_MANTIS",
  col = "lightgray"
)

boxplot(
  clin_final$MSI_SCORE_MANTIS,
  main = "Boxplot de MSI_SCORE_MANTIS",
  ylab = "MSI_SCORE_MANTIS",
  col= "lightgray"
)

# 4. Relação entre variáveis
#Idade vs. TMB
plot(clin_final$AGE, clin_final$TMB_NONSYNONYMOUS,
     col = "darkblue", pch = 16,
     main = "Idade vs TMB")

#Idade vs. MSI
plot(clin_final$AGE, clin_final$MSI_SCORE_MANTIS,
     col = "darkred", pch = 16)

#Stage vs. TMB
boxplot(TMB_NONSYNONYMOUS ~ AJCC_PATHOLOGIC_TUMOR_STAGE, data = clin_final,
        col = "lightgrey", las = 2)

#Sexo vs. TMB
boxplot(TMB_NONSYNONYMOUS ~ SEX, data = clin_final,
        col = c("pink", "lightblue"))

# 4. PCA da expressão
# Remover genes com variância zero
gene_var <- apply(expr_final, 1, var, na.rm = TRUE)
expr_pca <- expr_final[gene_var > 0, ]

# Usar apenas os genes mais variáveis (top 25%)
cutoff <- quantile(gene_var[gene_var > 0], 0.75)
expr_pca <- expr_pca[gene_var[gene_var > 0] >= cutoff, ]

cat("Dimensão da matriz usada para PCA:", dim(expr_pca), "\n")

# PCA com amostras nas linhas
pca <- prcomp(t(expr_pca), scale. = TRUE)

# Grupo de idade
clin_final$Grupo_Idade <- ifelse(clin_final$AGE < median(clin_final$AGE, na.rm = TRUE),
                                 "Jovem", "Idoso")

plot(
  pca$x[, 1], pca$x[, 2],
  col = as.factor(clin_final$Grupo_Idade),
  pch = 16,
  main = "PCA dos dados de expressão",
  xlab = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 1), "%)"),
  ylab = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 1), "%)")
)

legend(
  "bottomright",
  legend = levels(as.factor(clin_final$Grupo_Idade)),
  col = 1:length(levels(as.factor(clin_final$Grupo_Idade))),
  pch = 16,
  title = "Grupo Etário"
)