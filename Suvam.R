# ===========================================
# GSE42872 DE Analysis — Using fData(gset)
# ===========================================

rm(list=ls(all.names=TRUE)); gc()
if(!requireNamespace("BiocManager")) install.packages("BiocManager")
libs <- c("GEOquery","limma","pheatmap","ggplot2","dplyr")
for(p in libs){
  if(!requireNamespace(p)) BiocManager::install(p, ask=FALSE)
  library(p, character.only=TRUE)
}

# -----------------------------
# 1. Working directory
# -----------------------------
workdir <- "C:/Users/Suvam/Documents/GEO_analysis"
# <-- CHANGE to your folder
dir.create(workdir, showWarnings=FALSE, recursive=TRUE)
setwd(workdir)
dir.create("Results", showWarnings=FALSE)
cat("Working directory:", workdir, "\n")

# -----------------------------
# 2. Download GEO dataset
# -----------------------------
gset <- getGEO("GSE42872", GSEMatrix=TRUE, AnnotGPL=TRUE)[[1]]
expr_matrix <- exprs(gset)
pdata <- pData(gset)

# -----------------------------
# 3. Define groups dynamically
# -----------------------------
group <- ifelse(grepl("vemurafenib", pdata$source_name_ch1, ignore.case=TRUE), "treated",
                ifelse(grepl("vehicle|DMSO", pdata$source_name_ch1, ignore.case=TRUE), "control", NA))

keep <- !is.na(group)
expr_matrix <- expr_matrix[, keep, drop=FALSE]
group <- factor(group[keep])
pdata <- pdata[keep, ]

cat("Samples per group:\n")
print(table(group))

if(any(table(group) < 2)){
  stop("❌ One or more groups have fewer than 2 samples. Limma cannot run DE analysis.")
}

# -----------------------------
# 4. Map probes → gene symbols using fData(gset)
# -----------------------------
fdata <- fData(gset)

# Check for Gene Symbol column
symbol_col <- grep("Gene.*Symbol", colnames(fdata), ignore.case=TRUE, value=TRUE)[1]
if(is.na(symbol_col)) stop("❌ No Gene Symbol column found in fData(gset).")

# Create data frame with symbols
expr_df <- data.frame(expr_matrix)
expr_df$PROBEID <- rownames(expr_df)
expr_df$SYMBOL <- fdata[match(expr_df$PROBEID, rownames(fdata)), symbol_col]

# Remove probes with no symbol
expr_df <- expr_df[expr_df$SYMBOL!="" & !is.na(expr_df$SYMBOL), ]

# Average duplicates
expr_gene <- expr_df %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::summarise(dplyr::across(where(is.numeric), mean, na.rm=TRUE))

# -----------------------------
# 5. Convert to numeric matrix
# -----------------------------
expr_matrix2 <- as.matrix(expr_gene[, -1])
expr_matrix2 <- apply(expr_matrix2, 2, as.numeric)
rownames(expr_matrix2) <- expr_gene$SYMBOL
stopifnot(is.numeric(expr_matrix2))
cat("Expression matrix dimensions:", dim(expr_matrix2), "\n")

# -----------------------------
# 6. Limma DE analysis
# -----------------------------
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

fit <- lmFit(expr_matrix2, design)
contrast <- makeContrasts(treated_vs_control = treated - control, levels=design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

deg_all <- topTable(fit2, number=Inf, adjust.method="BH")
deg_up <- deg_all %>% dplyr::filter(logFC>1 & adj.P.Val<0.05)
deg_down <- deg_all %>% dplyr::filter(logFC< -1 & adj.P.Val<0.05)

# -----------------------------
# 7. Save DEG results
# -----------------------------
write.csv(deg_all, "Results/DEG_all.csv")
write.csv(deg_up, "Results/DEG_upregulated.csv")
write.csv(deg_down, "Results/DEG_downregulated.csv")

# -----------------------------
# 8. Volcano plot
# -----------------------------
deg_all$sign <- "NotSignificant"
deg_all$sign[deg_all$logFC>1 & deg_all$adj.P.Val<0.05] <- "Up"
deg_all$sign[deg_all$logFC< -1 & deg_all$adj.P.Val<0.05] <- "Down"

ggplot(deg_all, aes(logFC, -log10(adj.P.Val), color=sign)) +
  geom_point(alpha=0.7) +
  scale_color_manual(values=c("blue","grey","red")) +
  theme_minimal() +
  ggtitle("Volcano Plot: Treated vs Control") -> volcano
ggsave("Results/Volcano_plot.png", volcano, width=8, height=6)

# -----------------------------
# 9. Heatmap top 25 DEGs
# -----------------------------
top25 <- head(deg_all[order(deg_all$adj.P.Val),], 25)
top_expr <- expr_matrix2[rownames(expr_matrix2) %in% rownames(top25),]
annotation_col <- data.frame(Group=group)
rownames(annotation_col) <- colnames(top_expr)

pheatmap(top_expr, scale="row", annotation_col=annotation_col,
         filename="Results/Heatmap_top25_DEGs.png", main="Top 25 DEGs")

# -----------------------------
# 10. Summary
# -----------------------------
summary_text <- paste0(
  "Total DEGs: ", nrow(deg_all),
  "\nUpregulated: ", nrow(deg_up),
  "\nDownregulated: ", nrow(deg_down)
)
writeLines(summary_text, "Results/Summary.txt")

cat("✅ DE Analysis complete! Check the 'Results' folder.\n")
