# Registo de Alterações Fases 1 e 2 — Trabalho ECDB 2025/2026
## LUAD (TCGA PanCancer Atlas)

---

## 06_statistical_analysis.R

Volcano plot atualizado: genes upregulated a vermelho, downregulated a azul e não significativos a cinzento, com legenda. Anteriormente todos os genes significativos estavam a vermelho sem distinção.

Adicionadas contagens separadas de genes upregulated (logFC > 0) e downregulated (logFC < 0) entre estádios Early e Advanced.

Adicionada análise diferencial por sexo (modelo `~ SEX + AGE`). 

---

## 08_unsupervised_analysis.R

Top 500 genes para PCA/MDS/t-SNE/UMAP em vez dos 100 originais. O standard na literatura é 500 ou mais genes.

t-SNE (`Rtsne`) e UMAP (`uwot`) adicionados. Métodos ensinados na aula 5, ausentes no script original.

MDS clássico com `cmdscale` adicionado, com variância explicada por eixo calculada através dos eigenvalues.

Elbow method para k-means adicionado: loop k=2 a 10 com `sum(km$withinss)`, seguindo o código ensinado na aula 6. O script original não tinha justificação para a escolha de k=2.

`nstart=25` no k-means para garantir estabilidade do resultado.

`heatmap()` base substituído por `pheatmap` com anotações clínicas (Grupo_Idade, STAGE_GROUP, SEX).

Tabelas cruzadas dos clusters vs STAGE_GROUP e SEX adicionadas.

---

## fase1_report.Rmd

Análise diferencial por sexo adicionada entre o heatmap e o KEGG, com resultados e volcano plot.

Gráfico dos flat patterns adicionado na secção de pré-processamento, justificando visualmente o critério de filtro de variância pela mediana.

Correlações de Pearson e Spearman e gráfico `pairs()` adicionados na secção exploratória.

Boxplot de expressão por amostra adicionado como QC visual, complementar ao `plotDensities`.

---

## fase2_report.Rmd

MDS Leading LogFC, t-SNE e UMAP coloridos por sexo.

PCA colorida por sexo e por estádio adicionadas, além da versão por idade, permitindo comparar visualmente os três fatores.

Secção de k-means adicionada com elbow method e tabelas cruzadas.

