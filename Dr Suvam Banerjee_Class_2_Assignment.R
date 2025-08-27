classify_gene <- function(log2FC, padj) {
  if ((log2FC > 1) & (padj < 0.05)) {
    result <- "Upregulated"
  } else if ((log2FC < -1) & (padj < 0.05)) {
    result <- "Downregulated"
  } else {
    result <- "Not_Significant"
  }
  return(result)
}
input_dir <- "C:\\Users\\Suvam\\Downloads\\Raw_Data"
output_dir <- "Results"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

files_to_process <- c("DEGs_Data_1.csv", "DEGs_Data_2.csv")
deg_results <- list()

for (file in files_to_process) {
  cat("\nProcessing:", file, "\n")
  
  input_file <- file.path(input_dir, file)
  deg_data <- read.csv(input_file, header = TRUE)
  
  # replace NA padj with 1
  deg_data$padj[is.na(deg_data$padj)] <- 1
  
  # classify genes
  deg_data$status <- mapply(classify_gene, deg_data$logFC, deg_data$padj)
  
  # summary
  cat("Summary of gene status in", file, ":\n")
  print(table(deg_data$status))
  
  # save results
  deg_results[[file]] <- deg_data
  output_file <- file.path(output_dir, paste0("Processed_", file))
  write.csv(deg_data, output_file, row.names = FALSE)
  cat("Results saved to:", output_file, "\n")
}

