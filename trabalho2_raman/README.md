# Trabalho 2 — Extração de Conhecimento de Dados de Espectroscopia Raman

Análise de espectros Raman de compostos orgânicos puros, com o objetivo de explorar padrões espectrais, agrupar compostos por semelhança química e construir modelos de classificação automática. O trabalho foi desenvolvido em Python e organizado num notebook Jupyter.

## Compostos Analisados (Grupo 3)

| Composto | Família Química |
|---|---|
| Isobutylamine | Amina |
| Methanol | Álcool |
| Methyl isobutyl ketone | Cetona |
| N,N-Dimethylformamide | Amida |
| n-Heptane | Alcano linear |
| n-Hexane | Alcano linear |
| Propyl acetate | Éster |
| tert-Butanol | Álcool |
| tert-Butyl methyl ether | Éter |
| Tetrahydrofuran | Éter cíclico |
| Toluene | Aromático |

## Estrutura do Notebook

O notebook `trabalho2_grupo3.ipynb` está organizado nas seguintes secções:

1. **Carregamento e Exploração Inicial** — importação dos dados, identificação das variáveis de Raman Shift e estrutura do dataset
2. **Pré-processamento** — cropping espectral (300–2000 cm⁻¹), remoção de spikes (Whitaker-Hayes), suavização (Savitzky-Golay), correção de baseline (drPLS) e normalização Min-Max
3. **Estatística Descritiva e Exploração Gráfica** — espectros individuais por composto, identificação de bandas Raman relevantes e matriz de correlação espetral
4. **Análise Multivariada Não Supervisionada** — PCA, LDA e UMAP, com interpretação dos loadings e visualização dos agrupamentos
5. **Clustering** — Elbow Method para seleção de k, K-Means com k=11, validação com ARI e Silhouette Score, e matriz de contingência composto vs. cluster
6. **Aprendizagem Supervisionada** — comparação de Random Forest, SVM (RBF), KNN, Gradient Boosting e Regressão Logística com validação cruzada 5-fold; avaliação no conjunto de teste independente
7. **Importância das Variáveis** — Feature Importances do Random Forest, com interpretação das regiões espectrais mais discriminativas
8. **Discussão dos Resultados e Conclusão**

## Principais Resultados

- Pipeline de pré-processamento eficaz na eliminação de artefactos instrumentais e fluorescência
- Clustering K-Means com ARI = 0.856 e Silhouette = 0.761, reproduzindo com elevada fidelidade as classes reais sem supervisão
- Random Forest com acurácia de 100% no conjunto de teste, ancorada em bandas vibracionais quimicamente relevantes (região de fingerprint ~700–1100 cm⁻¹, respiração do anel aromático ~1000 cm⁻¹)

## Bibliotecas

O projeto foi desenvolvido em **Python 3**, com recurso às seguintes bibliotecas:

- `numpy`, `pandas` — manipulação de dados
- `scipy` — deteção de picos, suavização Savitzky-Golay
- `scikit-learn` — PCA, LDA, UMAP, K-Means, modelos de classificação, validação cruzada
- `matplotlib`, `seaborn` — visualização
- `pybaselines` — correção de baseline com drPLS

## Dados

Os dados de espectroscopia Raman não estão incluídos neste repositório. 

