# Get corresponding gene names from expression matrix (optimal solution for i = 3.3)
i <- 3.3
bladder_data <- loadRDS("data/bladder_cancer_clean.RDS")
expr_protein_coding <- bladder_data$X
filename <- paste0("data_res/sdib_km_res_bladder_clean_aggr_", i, ".RDS")
sdib_res <- readRDS(filename)
selected_indices <- which(sdib_res$weights > 0)
selected_genes <- colnames(expr_protein_coding)[selected_indices]

library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
ensembl_ids <- sub("\\..*", "", selected_genes)

gene_info <- getBM(
  attributes = c('ensembl_gene_id', 'external_gene_name', 'description'),
  filters = 'ensembl_gene_id',
  values = ensembl_ids,
  mart = ensembl
)

# Merge with weights
selected_df <- data.frame(
  ensembl_id = ensembl_ids,
  index = selected_indices,
  weight = sdib_res$weights[selected_indices]
)

selected_df <- merge(selected_df, gene_info, 
                     by.x = "ensembl_id", by.y = "ensembl_gene_id",
                     all.x = TRUE)

# Sort by weight (most important first)
selected_df <- selected_df[order(-selected_df$weight), ]

# More comprehensive marker lists based on Robertson (2017)
basal_markers <- c(
  "KRT5", "KRT6A", "KRT6B", "KRT6C", "KRT14", "KRT16", "KRT17",
  "EGFR", "CD44", "TP63", "DSC3", "S100A2", "S100P",
  "TBX2")

luminal_markers <- c(
  # Uroplakins
  "UPK1A", "UPK2", "UPK3A", "UPK3B",
  # Core transcription factors
  "GATA3", "FOXA1", "HNF1B", "GRHL3", "ELF3",
  # Signaling/receptors
  "PPARG", "FGFR3", "ERBB2", "ERBB3", "KRT20",
  # Other luminal-associated from Robertson (2017)
  "SNX31", "GPR160", "XBP1")

neuronal_markers <- c(
  "CHGA", "CHGB", "SYP", "NCAM1", "ENO2",
  "TUBB2B", "ELAVL4", "STMN2", "SNCG")

# Annotate
selected_df$subtype <- "Other"
selected_df$subtype[selected_df$external_gene_name %in% basal_markers] <- "Basal"
selected_df$subtype[selected_df$external_gene_name %in% luminal_markers] <- "Luminal"
selected_df$subtype[selected_df$external_gene_name %in% neuronal_markers] <- "Neuronal"

# Summary
print(table(selected_df$subtype))

# View all known markers
known_markers <- selected_df[selected_df$subtype != "Other", 
                             c("external_gene_name", "weight", "subtype")]
known_markers <- known_markers[order(known_markers$subtype, -known_markers$weight), ]
print(known_markers)

# By subtype
print(known_markers[known_markers$subtype == "Luminal", ])
print(known_markers[known_markers$subtype == "Neuronal", ])
print(known_markers[known_markers$subtype == "Basal", ])


library(ggplot2)
library(ggrepel)

selected_df$rank <- rank(-selected_df$weight)

ggplot(selected_df, aes(x = rank, y = weight, color = subtype)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(
    values = c("Basal" = "#E41A1C", "Luminal" = "#377EB8", 
               "Neuronal" = "#4DAF4A", "Other" = "gray80"),
    name = "Marker Type"
  ) +
  geom_text_repel(
    data = selected_df[selected_df$subtype != "Other", ],
    aes(label = external_gene_name),
    size = 3.5,
    fontface = "bold",
    max.overlaps = 20,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    segment.size = 0.2
  ) +
  scale_y_log10() +
  labs(
    title = "Sparse DIB Selected Genes Ranked by Weight",
    x = "Rank (1 = highest weight)",
    y = "Weight"
  ) +
  theme_bw() +
  theme(legend.position = "right")
