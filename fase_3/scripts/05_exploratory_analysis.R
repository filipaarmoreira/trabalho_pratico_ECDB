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


# Relações entre variáveis numéricas

boxplot(TMB_NONSYNONYMOUS ~ AJCC_PATHOLOGIC_TUMOR_STAGE, data = meta,
        col = "lightgrey", las = 2,
        main = "TMB por estádio")

boxplot(TMB_NONSYNONYMOUS ~ SEX, data = meta,
        col = c("pink", "lightblue"),
        main = "TMB por sexo")

# Correlações entre variáveis numéricas (Pearson e Spearman)

vars_cor <- c("AGE", "OS_MONTHS", "TMB_NONSYNONYMOUS",
              "ANEUPLOIDY_SCORE", "MSI_SCORE_MANTIS")
meta_num <- meta[, vars_cor]
meta_num <- meta_num[complete.cases(meta_num), ]

cat("\nCorrelações de Pearson entre variáveis clínicas:\n")
print(round(cor(meta_num, method = "pearson"), 3))

cat("\nCorrelações de Spearman entre variáveis clínicas:\n")
print(round(cor(meta_num, method = "spearman"), 3))

# Comparação visual de todos os pares de variáveis numéricas
pairs(
  meta_num,
  main = "Comparação de variáveis numéricas clínicas",
  pch  = 16,
  cex  = 0.4,
  col  = "steelblue"
)


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

# Boxplot de expressão por amostra (QC visual da normalização)
# usa subconjunto para legibilidade quando há muitas amostras
n_show   <- min(40, ncol(expr_mat))
idx_show <- round(seq(1, ncol(expr_mat), length.out = n_show))

boxplot(
  expr_mat[, idx_show],
  col     = "lightblue",
  las     = 2,
  outline = FALSE,
  main    = "Boxplot de expressão por amostra (subconjunto)",
  ylab    = "log2 expressão",
  xaxt    = "n"
)
abline(h = median(expr_mat, na.rm = TRUE), col = "red", lty = 2)

plotDensities(
  expr_mat,
  legend = FALSE,
  main = "Distribuição da expressão génica entre amostras"
)
