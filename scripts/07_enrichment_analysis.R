library(clusterProfiler)
library(pathview)
library(org.Hs.eg.db)
library(readr)
library(dplyr)

df <- read_csv("../data/processed/genes_significativos.csv")

# ver nomes das colunas
print(colnames(df))

# adaptar nomes se preciso
df2 <- df %>%
  select(Gene, logFC) %>%
  filter(!is.na(Gene), !is.na(logFC))

# converter IDs
conv <- bitr(df2$Gene,
             fromType = "SYMBOL",
             toType = "ENTREZID",
             OrgDb = org.Hs.eg.db)

df3 <- inner_join(conv, df2, by = c("SYMBOL" = "Gene"))

# criar vetor
gene_list <- df3$logFC
names(gene_list) <- df3$ENTREZID

# enrichment
kk <- enrichKEGG(gene = names(gene_list),
                 organism = "hsa")

# ver resultados
print(head(as.data.frame(kk)))

# plotar uma via
pathview(gene.data = gene_list,
         pathway.id = "04110",
         species = "hsa",
         limit = list(gene = 1))
