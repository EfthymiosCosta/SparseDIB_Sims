# Install BiocManager if needed
if (!require("BiocManager", quietly = TRUE))  install.packages("BiocManager")
#BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)
library(SummarizedExperiment)

# Load bladder dataset TCGA-BLCA
query <- GDCquery(
  project = "TCGA-BLCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# Download & prepare the data
GDCdownload(query)
data <- GDCprepare(query)
expr_matrix <- assay(data)
clinical_data <- colData(data)
expr_matrix_t <- t(expr_matrix)
dim(expr_matrix_t)  

# mRNA expression-based clusters - remove NA clusters
table(clinical_data$`paper_mRNA cluster`, useNA = "always")
molecular_subtypes <- clinical_data$`paper_mRNA cluster`
has_subtype <- !is.na(molecular_subtypes)

# Filter expression matrix & pre-processing
expr_matrix <- assay(data, "tpm_unstrand")
expr_matrix_t <- t(log2(expr_matrix + 1))
expr_clean <- expr_matrix_t[has_subtype, ]
subtypes_clean <- molecular_subtypes[has_subtype]
gene_info <- rowData(data)
colnames(gene_info)

# Check what gene types are available
table(gene_info$gene_type)
is_protein_coding <- gene_info$gene_type == "protein_coding"
expr_protein_coding <- expr_clean[, is_protein_coding]
dim(expr_protein_coding) 
gene_vars <- apply(expr_protein_coding, 2, var)

# A bunch of zero variance genes - need to remove those with variance e.g. < 0.01
housekeeping_genes <- which(gene_vars < 0.01)
expr_protein_coding <- expr_protein_coding[, -housekeeping_genes]
expr_protein_coding <- data.frame(expr_protein_coding)
dim(expr_protein_coding)
bladder_cancer_clean <- list(X = expr_protein_coding,
                             y = subtypes_clean)