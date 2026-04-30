# scripts/09_predictive_analysis.R
# Fase 3 – Análise Preditiva: Classificação de STAGE_GROUP (Early vs Advanced)
# Modelos: Elastic Net (glmnet), Random Forest, SVM
# Inclui: pré-processamento, validação cruzada, comparação de modelos,
#         seleção de atributos (filtro + embedded + RFE), importância de genes

library(Biobase)
library(limma)
library(caret)
library(glmnet)
library(randomForest)
library(e1071)
library(pROC)
library(ggplot2)
library(pheatmap)



# eset_ml é carregado pelo Rmd antes de correr este script

# ── 1. Preparar dados ──────────────────────────────────────────────────────────
set.seed(42)

# Target
y <- factor(pData(eset_ml)$STAGE_GROUP, levels = c("Early", "Advanced"))

# Split
train_idx <- createDataPartition(y, p = 0.75, list = FALSE)

y_train <- factor(pData(eset_ml)$STAGE_GROUP[train_idx], levels = c("Early", "Advanced"))
y_test  <- factor(pData(eset_ml)$STAGE_GROUP[-train_idx], levels = c("Early", "Advanced"))

# Feature selection (limma)
expr_train <- exprs(eset_ml)[, train_idx]
meta_train <- pData(eset_ml)[train_idx, ]

# Remover amostras com NA nas covariáveis
complete_idx <- complete.cases(meta_train[, c("AGE", "SEX")])
meta_train   <- meta_train[complete_idx, ]
expr_train   <- expr_train[, complete_idx]
y_train      <- y_train[complete_idx]

design <- model.matrix(~ STAGE_GROUP + AGE + SEX, data = meta_train)

fit <- lmFit(expr_train, design)
fit <- eBayes(fit)

results <- topTable(fit, coef = "STAGE_GROUPAdvanced", number = Inf)

top_genes <- rownames(results)[1:500]
top_genes <- top_genes[!is.na(top_genes)]

# Matrizes finais
X_train <- t(exprs(eset_ml)[top_genes, train_idx])
X_train <- X_train[complete_idx, ]
X_test  <- t(exprs(eset_ml)[top_genes, -train_idx])

# CV
ctrl <- trainControl(
  method          = "repeatedcv",
  number          = 5,
  repeats         = 3,
  classProbs      = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)
# ── 4. Modelo 1 – Elastic Net ────────────────────────────────────────────────

cat("\n=== Modelo 1: Elastic Net ===\n")

# Grid de hiperparâmetros 
tune_enet <- expand.grid(
  alpha  = c(0, 0.25, 0.5, 0.75, 1),
  lambda = 10^seq(-4, 0, length.out = 25)
)

set.seed(42)

model_glmnet <- train(
  x = X_train,
  y = y_train,
  method = "glmnet",
  metric = "ROC",
  trControl = ctrl,
  preProcess = c("center", "scale"),
  tuneGrid = tune_enet
)

# ── Melhores hiperparâmetros ────────────────────────────────────────────────

cat("Melhor alpha:", model_glmnet$bestTune$alpha,
    "| Melhor lambda:", round(model_glmnet$bestTune$lambda, 6), "\n")

# ── Genes selecionados (embedded feature selection) ─────────────────────────

coef_enet <- coef(model_glmnet$finalModel, s = model_glmnet$bestTune$lambda)

coef_df <- data.frame(
  Gene  = rownames(coef_enet),
  Coeff = as.numeric(coef_enet)
)

coef_df <- coef_df[
  coef_df$Gene != "(Intercept)" & coef_df$Coeff != 0,
]

coef_df <- coef_df[order(abs(coef_df$Coeff), decreasing = TRUE), ]

cat("Genes selecionados (coef != 0):", nrow(coef_df), "\n")

# ── Avaliação no test set ───────────────────────────────────────────────────

pred_enet <- predict(model_glmnet, X_test)

prob_enet <- predict(model_glmnet, X_test, type = "prob")[, "Advanced"]

cm_enet <- confusionMatrix(pred_enet, y_test, positive = "Advanced")

roc_enet <- roc(
  y_test,
  prob_enet,
  levels = c("Early", "Advanced"),
  direction = "<",
  quiet = TRUE
)

cat("Accuracy:", round(cm_enet$overall["Accuracy"], 3),
    "| AUC:", round(auc(roc_enet), 3), "\n")

# ── 5. Modelo 2 – Random Forest ─────────────────────────────────────────────

cat("\n=== Modelo 2: Random Forest ===\n")

# Grid de mtry (idealmente proporcional ao nº de features)
tune_rf <- expand.grid(
  mtry = unique(round(seq(5, ncol(X_train), length.out = 5)))
)

set.seed(42)

model_rf <- train(
  x          = X_train,
  y          = y_train,
  method     = "rf",
  trControl  = ctrl,
  tuneGrid   = tune_rf,
  metric     = "ROC",
  ntree      = 500,
  importance = TRUE
)

# Melhor hiperparâmetro
cat("Melhor mtry:", model_rf$bestTune$mtry, "\n")

# ── Avaliação no test set ───────────────────────────────────────────────

pred_rf <- predict(model_rf, X_test)

prob_rf <- predict(model_rf, X_test, type = "prob")[, "Advanced"]

cm_rf <- confusionMatrix(pred_rf, y_test, positive = "Advanced")

roc_rf <- roc(
  y_test,
  prob_rf,
  levels = c("Early", "Advanced"),
  direction = "<",
  quiet = TRUE
)

cat("Accuracy:", round(cm_rf$overall["Accuracy"], 3),
    "| AUC:", round(auc(roc_rf), 3), "\n")


# ── 6. Modelo 3 – SVM (RBF kernel) ───────────────────────────────────────

cat("\n=== Modelo 3: SVM (RBF kernel) ===\n")

# Grid mais robusto (melhor exploração)
tune_svm <- expand.grid(
  C     = 10^seq(-2, 2, length.out = 5),
  sigma = 10^seq(-4, -1, length.out = 5)
)

set.seed(42)

model_svm <- train(
  x          = X_train,
  y          = y_train,
  method     = "svmRadial",
  trControl  = ctrl,
  tuneGrid   = tune_svm,
  metric     = "ROC",
  preProcess = c("center", "scale")
)

# Melhor hiperparâmetros
cat("Melhor C:", model_svm$bestTune$C,
    "| Melhor sigma:", model_svm$bestTune$sigma, "\n")

# ── Avaliação no test set ───────────────────────────────────────────────

pred_svm <- predict(model_svm, X_test)

prob_svm <- predict(model_svm, X_test, type = "prob")[, "Advanced"]

cm_svm <- confusionMatrix(pred_svm, y_test, positive = "Advanced")

roc_svm <- roc(
  y_test,
  prob_svm,
  levels = c("Early", "Advanced"),
  direction = "<",
  quiet = TRUE
)

cat("Accuracy:", round(cm_svm$overall["Accuracy"], 3),
    "| AUC:", round(auc(roc_svm), 3), "\n")

cat("\n=== Comparação dos modelos (CV) ===\n")

resamps <- resamples(list(
  ElasticNet   = model_glmnet,
  RandomForest = model_rf,
  SVM          = model_svm
))

print(summary(resamps))

diffs <- diff(resamps)

cat("\nTeste estatístico de diferenças:\n")
print(summary(diffs))

bwplot(resamps, metric = "ROC",
       main = "Comparação de modelos – AUC (validação cruzada)")


# ── 8. Tabela resumo ──────────────────────────────────────────────────────────

perf_table <- data.frame(
  Modelo       = c("Elastic Net", "Random Forest", "SVM (RBF)"),
  Accuracy     = round(c(cm_enet$overall["Accuracy"],
                         cm_rf$overall["Accuracy"],
                         cm_svm$overall["Accuracy"]), 3),
  Sensitivity  = round(c(cm_enet$byClass["Sensitivity"],
                         cm_rf$byClass["Sensitivity"],
                         cm_svm$byClass["Sensitivity"]), 3),
  Specificity  = round(c(cm_enet$byClass["Specificity"],
                         cm_rf$byClass["Specificity"],
                         cm_svm$byClass["Specificity"]), 3),
  Balanced_Acc = round(c(cm_enet$byClass["Balanced Accuracy"],
                         cm_rf$byClass["Balanced Accuracy"],
                         cm_svm$byClass["Balanced Accuracy"]), 3),
  AUC          = round(c(auc(roc_enet), auc(roc_rf), auc(roc_svm)), 3)
)

rownames(perf_table) <- NULL

cat("\n=== Tabela resumo de desempenho ===\n")
print(perf_table)


# ── 9. Curvas ROC ─────────────────────────────────────────────────────────────

plot(roc_enet, col = "steelblue", lwd = 2,
     main = "Curvas ROC – conjunto de teste")

lines(roc_rf,  col = "firebrick", lwd = 2)
lines(roc_svm, col = "darkgreen", lwd = 2)

abline(a = 0, b = 1, lty = 2, col = "gray60")

legend("bottomright",
       legend = c(
         paste0("Elastic Net   (AUC = ", round(auc(roc_enet), 3), ")"),
         paste0("Random Forest (AUC = ", round(auc(roc_rf),  3), ")"),
         paste0("SVM (RBF)     (AUC = ", round(auc(roc_svm), 3), ")")
       ),
       col = c("steelblue","firebrick","darkgreen"),
       lwd = 2, bty = "n")


# ── 10. Importância de genes ──────────────────────────────────────────────────

## Random Forest
imp_rf <- varImp(model_rf, scale = TRUE)$importance
imp_rf$Overall <- rowMeans(imp_rf)
imp_rf$Gene <- rownames(imp_rf)
imp_rf <- imp_rf[order(imp_rf$Overall, decreasing = TRUE), ]

top20_rf <- head(imp_rf, 20)

cat("\n=== Top 20 genes – Random Forest ===\n")
print(top20_rf)

par(mar = c(4, 10, 3, 2))
barplot(rev(top20_rf$Overall), names.arg = rev(top20_rf$Gene),
        horiz = TRUE, las = 1, cex.names = 0.75, col = "firebrick",
        main = "Top 20 genes – Random Forest",
        xlab = "Importância")
par(mar = c(5, 4, 4, 2))


## Elastic Net
cat("\n=== Top 20 genes – Elastic Net ===\n")

top20_enet <- head(coef_df, 20)
print(top20_enet)

par(mar = c(4, 10, 3, 2))
barplot(rev(abs(top20_enet$Coeff)), names.arg = rev(top20_enet$Gene),
        horiz = TRUE, las = 1, cex.names = 0.75,
        col = ifelse(rev(top20_enet$Coeff) > 0, "steelblue", "darkorange"),
        main = "Top 20 genes – Elastic Net",
        xlab = "|Coeficiente|")

legend("bottomright",
       legend = c("Upregulated em Advanced","Downregulated em Advanced"),
       fill = c("steelblue","darkorange"), bty = "n", cex = 0.8)

par(mar = c(5, 4, 4, 2))


## Genes em consenso
top50_rf   <- head(imp_rf, 50)$Gene
top50_enet <- head(coef_df, 50)$Gene

genes_comuns <- intersect(top50_rf, top50_enet)

cat("\nGenes no top-50 dos dois modelos (", length(genes_comuns), "):\n")
print(genes_comuns)


# ── 11. RFE ──────────────────────────────────────────────────────────────────

cat("\n=== Seleção de atributos: RFE ===\n")

ctrl_rfe <- rfeControl(functions = rfFuncs, method = "cv", number = 5)

top200 <- head(top_genes, 200)

X_rfe <- X_train[, colnames(X_train) %in% top200]

rfe_result <- rfe(
  x          = X_rfe,
  y          = y_train,
  sizes      = c(10, 20, 30, 50, 100),
  rfeControl = ctrl_rfe
)

cat("Número ótimo de genes:", rfe_result$optsize, "\n")
cat("Accuracy CV:", round(max(rfe_result$results$Accuracy), 3), "\n")

genes_rfe <- predictors(rfe_result)
print(head(genes_rfe, 30))

plot(rfe_result, type = c("g", "o"),
     main = "RFE – Performance vs número de genes")


# ── 12. Heatmap ──────────────────────────────────────────────────────────────

genes_heat <- if (length(genes_comuns) >= 5) {
  genes_comuns[1:min(30, length(genes_comuns))]
} else {
  top50_rf[1:min(30, length(top50_rf))]
}

genes_heat <- genes_heat[genes_heat %in% rownames(eset_ml)]

annotation_col <- data.frame(
  STAGE_GROUP = pData(eset_ml)$STAGE_GROUP,
  SEX         = pData(eset_ml)$SEX
)

rownames(annotation_col) <- colnames(eset_ml)

ann_colors <- list(
  STAGE_GROUP = c(Early = "steelblue", Advanced = "firebrick"),
  SEX         = c(Male = "lightblue", Female = "pink")
)

pheatmap(
  exprs(eset_ml)[genes_heat, ],
  scale = "row",
  show_colnames = FALSE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  main = "Heatmap – genes mais importantes"
)


# ── 13. Guardar ──────────────────────────────────────────────────────────────

dir.create("../data/processed", recursive = TRUE, showWarnings = FALSE)

write.csv(perf_table, "../data/processed/model_performance.csv", row.names = FALSE)

saveRDS(model_glmnet, "../data/processed/model_elasticnet.rds")
saveRDS(model_rf,     "../data/processed/model_rf.rds")
saveRDS(model_svm,    "../data/processed/model_svm.rds")
saveRDS(rfe_result,   "../data/processed/rfe_result.rds")
saveRDS(perf_table,   "../data/processed/perf_table.rds")

write.csv(imp_rf,  "../data/processed/importance_rf.csv",   row.names = FALSE)
write.csv(coef_df, "../data/processed/importance_enet.csv", row.names = FALSE)

write.csv(
  data.frame(
    Gene            = genes_comuns,
    RF_rank         = match(genes_comuns, imp_rf$Gene),
    ElasticNet_rank = match(genes_comuns, coef_df$Gene)
  ),
  "../data/processed/genes_consenso.csv",
  row.names = FALSE
)
