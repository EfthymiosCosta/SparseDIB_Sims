library(VarSelLCM)

# Create output directory
if (!dir.exists("data_res")) dir.create("data_res")

bladder_cancer <- readRDS('data/bladder_cancer_clean.rds')
labels <- bladder_cancer$y
X <- as.data.frame(bladder_cancer$X)
dim(X)

labels_aggr <- as.character(labels)
labels_aggr[grepl("Luminal", labels_aggr)] <- "Luminal"

K <- length(unique(labels_aggr))

cat("Running VarSelCluster with K =", K, "clusters...\n")

vslcm <- VarSelCluster(x = X, gvals = K,
                       initModel = 1000)

result <- list(
  clustering = vslcm@partitions@zMAP,
  relevant_features = vslcm@model@names.relevant
)

output_file <- "data_res/varselcluster_bladder.RDS"
saveRDS(result, output_file)

cat("VarSelCluster analysis complete!\n")
cat("Clustering saved to:", output_file, "\n")
cat("Number of relevant features:", length(result$relevant_features), "\n")