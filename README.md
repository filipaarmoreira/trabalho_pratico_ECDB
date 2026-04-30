# Análise de Dados Ómicos: Adenocarcinoma Pulmonar (LUAD)
Este projeto consiste numa análise bioinformática de dados de expressão génica de Adenocarcinoma Pulmonar (LUAD), utilizando o estudo Lung Adenocarcinoma (TCGA PanCancer Atlas), obtido através do **cBioPortal**.
O trabalho foi desenvolvido em **R** e está dividido em três fases principais: preparação e exploração dos dados, análise não supervisionada e análise preditiva.

## Autores
-Ana Filipa Moreira

-Diana Oliveira

-Márcia Silva

-Rodolfo Ferreira

## Estrutura do Projeto
O trabalho está organizado em três fases sequenciais:

### Fase 1: Preparação, Exploração e Análise Diferencial
**Objetivo:** preparar, integrar e explorar os dados clínicos e de expressão génica, realizando também análise diferencial entre grupos tumorais.

**Principais Tarefas:**

- Importação dos dados clínicos e de expressão génica;

- Pré-processamento e integração dos dados num `ExpressionSet`;

- Análise exploratória inicial;

- Análise estatística univariada;

- Análise de expressão diferencial entre tumores Early e Advanced;

- Análise de enriquecimento funcional.

### Fase 2: Análise Não Supervisionada
**Objetivo:** Identificar padrões e agrupamentos naturais nos dados sem usar etiquetas prévias.

**Principais tarefas:**

- Redução de dimensionalidade com PCA, MDS, t-SNE e UMAP;

- Clustering hierárquico e k-means;

- Visualização com heatmaps e dendrogramas;

- Avaliação da associação dos agrupamentos com variáveis clínicas, como idade, sexo e estádio tumoral.

### Fase 3: Análise Preditiva
**Objetivo:** Construir e avaliar modelos de Machine Learning para classificação e previsão baseada nos dados ómicos.

**Principais tarefas:**

- Treino e teste de modelos preditivos;

- Aplicação de Elastic Net, Random Forest e Support Vector Machine;

- Seleção de genes relevantes através de limma, Elastic Net e RFE;

- Avaliação da performance com métricas como matriz de confusão, sensibilidade, especificidade, accuracy e AUC;

- Identificação de potenciais biomarcadores preditivos.

## Tecnologias e Bibliotecas Utilizadas
O projeto foi desenvolvido em **R**, com recurso a bibliotecas especializadas para análise estatística, visualização e machine learning.

Principais ferramentas utilizadas:

- cBioPortal — obtenção dos dados TCGA PanCancer Atlas;
- Biobase — criação e manipulação do `ExpressionSet`;
- limma — análise de expressão diferencial;
- ggplot2 — visualização de dados;
- pheatmap / ComplexHeatmap — construção de heatmaps;
- caret — treino e avaliação de modelos preditivos;
- randomForest / glmnet / e1071 — modelos de machine learning.

## Como Executar

1. Certificar que o R e o RStudio estão instalados.
2. Instalar as dependências indicadas nos ficheiros `.R` ou `.Rmd`.
3. Executar os scripts pela ordem indicada em cada fase.
4. Abrir os relatórios HTML para consultar os resultados detalhados.

