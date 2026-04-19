# scripts/05_exploratory_analysis.R

library(Biobase)
library(limma)

# Extrair dados do ExpressionSet

expr_mat <- exprs(eset)
meta <- pData(eset)

cat("Dimensão do dataset clínico final:", dim(meta), "\n")
cat("Dimensão da matriz de expressão final:", dim(expr_mat), "\n")

# Variáveis categóricas clínicas

tabela_estadios <- table(meta$AJCC_PATHOLOGIC_TUMOR_STAGE, useNA = "ifany")
print(tabela_estadios)

barplot(
  tabela_estadios,
  main = "Distribuição dos estádios do cancro",
  col = "coral",
  las = 2,
  cex.names = 0.7
)

tabela_sexo <- table(meta$SEX, useNA = "ifany")
print(tabela_sexo)

barplot(
  tabela_sexo,
  main = "Distribuição por sexo",
  col = "lightpink"
)

tabela_os <- table(meta$OS_STATUS, useNA = "ifany")
print(tabela_os)

barplot(
  tabela_os,
  main = "Distribuição de OS_STATUS",
  col = "lightgreen"
)

tabela_pfs <- table(meta$PFS_STATUS, useNA = "ifany")
print(tabela_pfs)

barplot(
  tabela_pfs,
  main = "Distribuição de PFS_STATUS",
  col = "lightblue"
)


# Variáveis numéricas clínicas e moleculares

summary(meta$AGE)
hist(meta$AGE, main = "Distribuição da idade dos pacientes",
     xlab = "Idade (anos)", ylab = "Frequência", col = "lightblue")
boxplot(meta$AGE, main = "Boxplot da idade", ylab = "Idade", col = "lightblue")

summary(meta$OS_MONTHS)
hist(meta$OS_MONTHS, main = "Distribuição de OS_MONTHS",
     xlab = "Overall survival (meses)", col = "lightgreen")
boxplot(meta$OS_MONTHS, main = "Boxplot de OS_MONTHS",
        ylab = "Meses", col = "lightgreen")

summary(meta$ANEUPLOIDY_SCORE)
hist(meta$ANEUPLOIDY_SCORE, main = "Distribuição de Aneuploidy Score",
     xlab = "Score", col = "lightyellow")
boxplot(meta$ANEUPLOIDY_SCORE, main = "Boxplot - Aneuploidy Score",
        ylab = "Score", col = "lightyellow")

summary(meta$TMB_NONSYNONYMOUS)
hist(meta$TMB_NONSYNONYMOUS, main = "Distribuição de TMB",
     xlab = "TMB", col = "lavender")
boxplot(meta$TMB_NONSYNONYMOUS, main = "Boxplot de TMB",
        ylab = "TMB", col = "lavender")

summary(meta$MSI_SCORE_MANTIS)
hist(meta$MSI_SCORE_MANTIS, main = "Distribuição de MSI_SCORE_MANTIS",
     xlab = "MSI_SCORE_MANTIS", col = "lightgray")
boxplot(meta$MSI_SCORE_MANTIS, main = "Boxplot de MSI_SCORE_MANTIS",
        ylab = "MSI_SCORE_MANTIS", col = "lightgray")


# Relações simples entre variáveis

plot(meta$AGE, meta$TMB_NONSYNONYMOUS,
     col = "darkblue", pch = 16,
     main = "Idade vs TMB",
     xlab = "Idade", ylab = "TMB")

plot(meta$AGE, meta$MSI_SCORE_MANTIS,
     col = "darkred", pch = 16,
     main = "Idade vs MSI",
     xlab = "Idade", ylab = "MSI_SCORE_MANTIS")

boxplot(TMB_NONSYNONYMOUS ~ AJCC_PATHOLOGIC_TUMOR_STAGE, data = meta,
        col = "lightgrey", las = 2,
        main = "TMB por estádio")

boxplot(TMB_NONSYNONYMOUS ~ SEX, data = meta,
        col = c("pink", "lightblue"),
        main = "TMB por sexo")


# QC da expressão

sample_means <- colMeans(expr_mat, na.rm = TRUE)
sample_medians <- apply(expr_mat, 2, median, na.rm = TRUE)
sample_sds <- apply(expr_mat, 2, sd, na.rm = TRUE)

cat("\nResumo das medianas por amostra:\n")
print(summary(sample_medians))

cat("\nResumo das médias por amostra:\n")
print(summary(sample_means))

cat("\nResumo dos desvios padrão por amostra:\n")
print(summary(sample_sds))

plotDensities(
  expr_mat,
  legend = FALSE,
  main = "Distribuição da expressão génica entre amostras"
)
