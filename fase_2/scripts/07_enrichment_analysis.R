# scripts/07_enrichment_analysis.R

library(clusterProfiler)
library(pathview)
library(org.Hs.eg.db)
library(readr)
library(dplyr)

# Pastas de output

dir.create("../results", showWarnings = FALSE, recursive = TRUE)
dir.create("../data/processed", showWarnings = FALSE, recursive = TRUE)

# Ler resultados do limma

res <- read_csv(
  "../data/processed/resultados_limma_stage_age_sex.csv",
  show_col_types = FALSE
)

stopifnot(all(c("Gene", "logFC", "P.Value", "adj.P.Val") %in% colnames(res)))

res <- res %>%
  select(Gene, logFC, P.Value, adj.P.Val) %>%
  filter(!is.na(Gene), !is.na(logFC)) %>%
  distinct(Gene, .keep_all = TRUE)


# Preparar IDs dos genes
#    - genes "normais" vão por SYMBOL
#    - genes "ENTREZ_xxx" usam fallback

res <- res %>%
  mutate(
    GENE_SYMBOL = ifelse(grepl("^ENTREZ_", Gene), NA_character_, Gene),
    ENTREZ_FALLBACK = ifelse(
      grepl("^ENTREZ_", Gene),
      sub("^ENTREZ_", "", Gene),
      NA_character_
    )
  )

symbols_to_map <- unique(na.omit(res$GENE_SYMBOL))

conv <- bitr(
  symbols_to_map,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

res2 <- res %>%
  left_join(conv, by = c("GENE_SYMBOL" = "SYMBOL")) %>%
  mutate(
    ENTREZID = ifelse(is.na(ENTREZID), ENTREZ_FALLBACK, ENTREZID),
    ENTREZID = as.character(ENTREZID),
    logFC = as.numeric(logFC)
  ) %>%
  filter(!is.na(ENTREZID), ENTREZID != "", !is.na(logFC))

cat("Genes com ENTREZ após conversão:", nrow(res2), "\n")

# Universo e genes significativos

universe_entrez <- unique(res2$ENTREZID)

sig_entrez <- res2 %>%
  filter(adj.P.Val < 0.05) %>%
  pull(ENTREZID) %>%
  unique()

cat("Genes no universo:", length(universe_entrez), "\n")
cat("Genes significativos com ENTREZ:", length(sig_entrez), "\n")

if (length(sig_entrez) == 0) {
  stop("Não há genes significativos com ENTREZID para correr enrichKEGG().")
}


# Vetor de logFC para o pathviews

gene_df <- res2 %>%
  group_by(ENTREZID) %>%
  summarise(logFC = mean(logFC, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(logFC), !is.na(ENTREZID), ENTREZID != "")

gene_list <- gene_df$logFC
names(gene_list) <- gene_df$ENTREZID

cat("gene_list é numérico:", is.numeric(gene_list), "\n")
cat("NAs em gene_list:", sum(is.na(gene_list)), "\n")
cat("NAs nos nomes:", sum(is.na(names(gene_list))), "\n")
cat("Nomes vazios:", sum(names(gene_list) == ""), "\n")
cat("Duplicados nos nomes:", anyDuplicated(names(gene_list)), "\n")


# Enrichment KEGG

kk <- enrichKEGG(
  gene = sig_entrez,
  universe = universe_entrez,
  organism = "hsa",
  pAdjustMethod = "BH"
)

kk_df <- as.data.frame(kk)

write_csv(
  kk_df,
  "../data/processed/kegg_enrichment_stage_age_sex.csv"
)

cat("\nTop vias KEGG:\n")
print(head(kk_df))

if (nrow(kk_df) > 0) {
  pdf("../results/kegg_dotplot_stage_age_sex.pdf", width = 9, height = 6)
  print(dotplot(kk, showCategory = min(15, nrow(kk_df))))
  dev.off()
}

# Gerar a via hsa04110 (Cell cycle)

chosen_id <- "04110"
chosen_desc <- "Cell cycle"

cat("\nVia escolhida para pathview:", chosen_id, "-", chosen_desc, "\n")

# tentar versão nativa primeiro
pv_out <- tryCatch(
  {
    pathview(
      gene.data = gene_list,
      pathway.id = chosen_id,
      species = "hsa",
      gene.idtype = "entrez",
      kegg.native = TRUE,
      out.suffix = "stage_age_sex"
    )
  },
  error = function(e) {
    cat("Falhou kegg.native=TRUE. Tentar kegg.native=FALSE...\n")
    pathview(
      gene.data = gene_list,
      pathway.id = chosen_id,
      species = "hsa",
      gene.idtype = "entrez",
      kegg.native = FALSE,
      out.suffix = "stage_age_sex"
    )
  }
)

cat("Mapa KEGG gerado para a via:", chosen_id, "-", chosen_desc, "\n")

