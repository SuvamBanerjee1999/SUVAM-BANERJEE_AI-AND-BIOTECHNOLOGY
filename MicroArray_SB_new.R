# =====================================================================
#   Microarray Preprocessing Workflow for GSE79973
# =====================================================================

# 0. Install & load packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", quiet = TRUE)
BiocManager::install(c("GEOquery","affy","arrayQualityMetrics","genefilter"), ask = FALSE, update = TRUE)
install.packages(c("dplyr","matrixStats"), dependencies = TRUE, quiet = TRUE)

library(GEOquery)
library(affy)
library(arrayQualityMetrics)
library(dplyr)
library(matrixStats)
library(genefilter)

# 1. Load series matrix (use local file first)
local_file <- "GSE79973_series_matrix.txt.gz"
if (file.exists(local_file)) {
  cat("Loading local series matrix...\n")
  gse_data <- getGEO(filename = local_file)
} else {
  cat("Local series matrix not found, downloading safely...\n")
  gse_data <- getGEO("GSE79973", GSEMatrix = TRUE, getGPL = FALSE)
}

expression_data <- exprs(gse_data[[1]])
feature_data    <- fData(gse_data[[1]])
phenotype_data  <- pData(gse_data[[1]])

# 2. Optional: Read CEL files if present
cel_dir <- "Raw_Data/CEL_Files"
cel_files <- list.files(cel_dir, pattern = "\\.CEL", full.names = TRUE, ignore.case = TRUE)

if (length(cel_files) > 0) {
  cat("CEL files found. Reading raw data...\n")
  raw_data <- ReadAffy(celfile.path = cel_dir)
  
  # QC before normalization
  arrayQualityMetrics(raw_data, outdir = "Results/QC_Raw_Data", force = TRUE, do.logtransform = TRUE)
  
  # RMA normalization
  normalized_data <- rma(raw_data)
  
  # QC after normalization
  arrayQualityMetrics(normalized_data, outdir = "Results/QC_Normalized_Data", force = TRUE)
  
  processed_data <- as.data.frame(exprs(normalized_data))
  
} else {
  cat("No CEL files found. Using series matrix expression data only.\n")
  processed_data <- as.data.frame(expression_data)
}

# 3. Filter low-intensity probes
row_median <- rowMedians(as.matrix(processed_data))
threshold <- 3.5
processed_data <- processed_data[row_median > threshold, ]
colnames(processed_data) <- rownames(phenotype_data)

# 4. Define groups (normal vs cancer)
groups <- factor(
  phenotype_data$source_name_ch1,
  levels = c("gastric mucosa","gastric adenocarcinoma"),
  labels = c("normal","cancer")
)
phenotype_data$group <- groups

# 5. Save processed data
dir.create("Results", showWarnings = FALSE)
write.csv(processed_data, "Results/Processed_Expression_Data.csv", row.names = TRUE)
write.csv(phenotype_data, "Results/Phenotype_Data.csv", row.names = TRUE)

# 6. Summary output
cat("\n===== Workflow Summary =====\n")
cat("Transcripts after filtering:", nrow(processed_data), "\n")
cat("Sample groups:\n")
print(table(phenotype_data$group))
cat("\n✅ Workflow completed successfully.\n")
pdf("Boxplot_Normalized.pdf")
boxplot(processed_data, main="Boxplot After Normalization", las=2, col="lightblue")
dev.off()
pca_res <- prcomp(t(processed_data), scale.=TRUE)
pdf("PCA_Normalized.pdf")
plot(pca_res$x[,1:2], col=as.numeric(phenotype_data$group), pch=19,
     xlab="PC1", ylab="PC2", main="PCA After Normalization")
legend("topright", legend=levels(phenotype_data$group), col=1:2, pch=19)
dev.off()
nrow(expression_data)
