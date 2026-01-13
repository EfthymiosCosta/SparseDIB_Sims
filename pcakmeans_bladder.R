# Create output directory
if (!dir.exists("data_res")) dir.create("data_res")

bladder_cancer <- readRDS('data/bladder_cancer_clean.rds')
labels <- bladder_cancer$y
X <- as.data.frame(bladder_cancer$X)
dim(X)

labels_aggr <- as.character(labels)
labels_aggr[grepl("Luminal", labels_aggr)] <- "Luminal"
K <- length(unique(labels_aggr))

cat("Running PCA + K-Means with K =", K, "clusters...\n")

pca_result <- prcomp(X, center = TRUE, scale. = TRUE)

cat("Clustering on first PC...\n")
pc1_projection <- pca_result$x[, 1, drop = FALSE]
kmeans_pc1 <- kmeans(pc1_projection, centers = K, nstart = 1000)
output_file_pc1 <- "data_res/pca1_kmeans_bladder.RDS"
saveRDS(kmeans_pc1$cluster, output_file_pc1)

cat("Clustering on first 2 PCs...\n")
pc2_projection <- pca_result$x[, 1:2]
kmeans_pc2 <- kmeans(pc2_projection, centers = K, nstart = 1000)
output_file_pc2 <- "data_res/pca2_kmeans_bladder.RDS"
saveRDS(kmeans_pc2$cluster, output_file_pc2)

cat("\nPCA K-means analysis complete!\n")