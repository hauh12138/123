rm(list = ls())
gc()
options(stringsAsFactors = FALSE, timeout = 1200)
setwd("E:\\HuaweiMoveData\\Users\\潘懿\\Desktop\\abab")
############################################################
## 0. Packages and directories
############################################################

need_pkgs <- c(
  "TCGAbiolinks", "SummarizedExperiment", "DESeq2", "WGCNA",
  "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi",
  "ggplot2", "pheatmap", "data.table", "dplyr", "tibble",
  "survival", "survminer", "pROC", "Seurat", "GEOquery", "Matrix", "tidyr"
)

missing_pkgs <- need_pkgs[!vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "),
       "\nPlease install them before running this pipeline.")
}

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(DESeq2)
  library(WGCNA)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
  library(pheatmap)
  library(data.table)
  library(dplyr)
  library(tibble)
  library(survival)
  library(survminer)
  library(pROC)
  library(Seurat)
  library(GEOquery)
  library(Matrix)
})

allowWGCNAThreads()

make_dirs <- function() {
  dirs <- c(
    "data", "data/raw", "data/processed", "data/clinical", "data/scRNA",
    "results", "results/tables", "results/figures",
    "results/figures/DEG", "results/figures/WGCNA", "results/figures/enrichment",
    "results/figures/survival", "results/figures/scRNA",
    "results/RData"
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}
make_dirs()

clean_ensembl <- function(x) sub("\\..*$", "", x)
`%||%` <- function(a, b) if (!is.null(a)) a else b

save_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, file, row.names = FALSE)
}

safe_ggsave <- function(filename, plot, width = 7, height = 6, dpi = 600) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename, plot = plot, width = width, height = height, dpi = dpi)
}

collapse_by_symbol <- function(mat, anno, stat = c("mean", "sum")) {
  stat <- match.arg(stat)
  anno2 <- anno %>%
    dplyr::select(ensembl_clean, gene_name) %>%
    dplyr::filter(!is.na(gene_name), gene_name != "") %>%
    dplyr::distinct()
  row_df <- data.frame(ensembl_clean = clean_ensembl(rownames(mat)), row_index = seq_len(nrow(mat)))
  map <- dplyr::left_join(row_df, anno2, by = "ensembl_clean") %>%
    dplyr::filter(!is.na(gene_name), gene_name != "")
  mat2 <- mat[map$row_index, , drop = FALSE]
  rownames(mat2) <- map$gene_name
  gene_groups <- split(seq_len(nrow(mat2)), rownames(mat2))
  out <- lapply(gene_groups, function(idx) {
    if (length(idx) == 1) return(mat2[idx, ])
    if (stat == "sum") colSums(mat2[idx, , drop = FALSE]) else colMeans(mat2[idx, , drop = FALSE])
  })
  out <- do.call(rbind, out)
  rownames(out) <- names(gene_groups)
  out
}

############################################################
## 1. TCGA-COAD bulk RNA-seq download / loading
############################################################

message("========== Step 1: TCGA-COAD bulk RNA-seq ==========")
se_file <- "data/raw/TCGA_COAD_SE_STAR_Counts.rds"

if (!file.exists(se_file)) {
  query <- GDCquery(
    project = "TCGA-COAD",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  GDCdownload(query, method = "api", files.per.chunk = 20)
  se <- GDCprepare(query)
  saveRDS(se, se_file)
} else {
  se <- readRDS(se_file)
}

assay_name <- if ("unstranded" %in% assayNames(se)) "unstranded" else assayNames(se)[1]
count_matrix <- assay(se, assay_name)
clinical <- as.data.frame(colData(se))
gene_anno <- as.data.frame(rowData(se))

gene_anno$ensembl_clean <- clean_ensembl(gene_anno$gene_id %||% rownames(gene_anno))
if (!"gene_name" %in% colnames(gene_anno)) gene_anno$gene_name <- rownames(gene_anno)
if (!"gene_type" %in% colnames(gene_anno)) gene_anno$gene_type <- NA_character_

clinical$sample_barcode <- colnames(count_matrix)
clinical$patient_barcode <- substr(clinical$sample_barcode, 1, 12)

# Prefer colData shortLetterCode if present. Fall back to TCGA barcode sample-type code.
if ("shortLetterCode" %in% colnames(clinical)) {
  sample_type <- clinical$shortLetterCode
} else {
  sample_type_code <- substr(colnames(count_matrix), 14, 15)
  sample_type <- ifelse(sample_type_code == "11", "NT", ifelse(sample_type_code == "01", "TP", NA))
}

keep_samples <- sample_type %in% c("TP", "NT")
count_matrix <- count_matrix[, keep_samples, drop = FALSE]
clinical <- clinical[keep_samples, , drop = FALSE]
clinical$shortLetterCode_final <- sample_type[keep_samples]
clinical$group <- factor(ifelse(clinical$shortLetterCode_final == "TP", "Tumor", "Normal"),
                         levels = c("Normal", "Tumor"))
rownames(clinical) <- colnames(count_matrix)

saveRDS(count_matrix, "data/processed/count_matrix_TP_NT.rds")
saveRDS(clinical, "data/clinical/clinical_TP_NT.rds")
saveRDS(gene_anno, "data/processed/gene_annotation.rds")

sample_summary <- clinical %>% dplyr::count(group, name = "n")
save_csv(sample_summary, "results/tables/Sample_Summary.csv")

############################################################
## 2. DESeq2 differential expression analysis
############################################################

message("========== Step 2: DESeq2 differential expression ==========")

before_filter_genes <- nrow(count_matrix)
keep_genes <- rowSums(count_matrix >= 10) >= max(3, ceiling(0.1 * ncol(count_matrix)))
count_filtered <- count_matrix[keep_genes, , drop = FALSE]
after_filter_genes <- nrow(count_filtered)

save_csv(data.frame(
  genes_before_filter = before_filter_genes,
  genes_after_filter = after_filter_genes,
  filter_rule = "count >= 10 in at least max(3, 10% samples)"
), "results/tables/Gene_Filter_Summary.csv")

coldata <- data.frame(
  row.names = colnames(count_filtered),
  group = clinical[colnames(count_filtered), "group"]
)

stopifnot(identical(colnames(count_filtered), rownames(coldata)))

dds <- DESeqDataSetFromMatrix(
  countData = round(count_filtered),
  colData = coldata,
  design = ~ group
)

dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", "Tumor", "Normal"))
res_df <- as.data.frame(res)
res_df$GeneID <- rownames(res_df)
res_df$ensembl_clean <- clean_ensembl(res_df$GeneID)
res_df <- dplyr::left_join(
  res_df,
  gene_anno %>% dplyr::select(ensembl_clean, gene_name, gene_type) %>% dplyr::distinct(),
  by = "ensembl_clean"
)
res_df <- res_df[!is.na(res_df$padj), ]
res_df$DEG_strict <- with(res_df, padj < 0.01 & abs(log2FoldChange) > 1.5)
res_df$DEG_loose <- with(res_df, padj < 0.05 & abs(log2FoldChange) > 1)
res_df$DEG_Group <- "Not significant"
res_df$DEG_Group[res_df$padj < 0.01 & res_df$log2FoldChange > 1.5] <- "Up"
res_df$DEG_Group[res_df$padj < 0.01 & res_df$log2FoldChange < -1.5] <- "Down"

save_csv(res_df, "results/tables/DEG_all_annotated.csv")
save_csv(res_df[res_df$DEG_strict, ], "results/tables/DEG_strict_annotated.csv")
save_csv(res_df[res_df$DEG_loose, ], "results/tables/DEG_loose_annotated.csv")

DEG_summary <- data.frame(
  total_tested_genes = nrow(res_df),
  strict_up = sum(res_df$padj < 0.01 & res_df$log2FoldChange > 1.5),
  strict_down = sum(res_df$padj < 0.01 & res_df$log2FoldChange < -1.5),
  loose_up = sum(res_df$padj < 0.05 & res_df$log2FoldChange > 1),
  loose_down = sum(res_df$padj < 0.05 & res_df$log2FoldChange < -1),
  strict_threshold = "padj < 0.01 and |log2FC| > 1.5",
  loose_threshold = "padj < 0.05 and |log2FC| > 1",
  method = "DESeq2 on TCGA-COAD STAR-Counts"
)
save_csv(DEG_summary, "results/tables/DEG_Summary.csv")

saveRDS(dds, "data/processed/dds_TP_NT.rds")
vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)
saveRDS(vsd, "data/processed/vsd_TP_NT.rds")
saveRDS(vsd_mat, "data/processed/vsd_matrix_TP_NT.rds")

############################################################
## Sample Pearson correlation heatmap
############################################################

library(pheatmap)
library(RColorBrewer)
vsd <- vst(dds)
mat <- assay(vsd)
dir.create("results/figures/QC", recursive = TRUE, showWarnings = FALSE)

# mat: genes x samples, from vst(dds)
# 为避免低变异基因稀释相关性，可选取变异度最高的前 5000 个基因
gene_var <- apply(mat, 1, var, na.rm = TRUE)

top_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(5000, length(gene_var))]

mat_top <- mat[top_genes, ]

# 计算样本间 Pearson 相关性
sample_cor <- cor(
  mat_top,
  method = "pearson",
  use = "pairwise.complete.obs"
)

# 样本分组注释
annotation_col <- data.frame(
  Group = coldata[colnames(sample_cor), "group"]
)

rownames(annotation_col) <- colnames(sample_cor)

# 保证行列注释一致
annotation_row <- annotation_col

# 保存相关性矩阵
write.csv(
  sample_cor,
  "results/tables/Sample_Pearson_Correlation.csv"
)

# PNG
png(
  filename = "results/figures/QC/Sample_Pearson_Correlation_Heatmap.png",
  width = 2400,
  height = 2200,
  res = 300
)

pheatmap(
  sample_cor,
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  show_colnames = FALSE,
  show_rownames = FALSE,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "average",
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  main = "Sample Pearson correlation heatmap",
  border_color = NA,
  fontsize = 10
)

dev.off()

# PDF
pdf(
  file = "results/figures/QC/Sample_Pearson_Correlation_Heatmap.pdf",
  width = 8,
  height = 7
)

pheatmap(
  sample_cor,
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  show_colnames = FALSE,
  show_rownames = FALSE,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "average",
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  main = "Sample Pearson correlation heatmap",
  border_color = NA,
  fontsize = 10
)

dev.off()

# PCA
# 补充：prcomp 结果需手动计算方差解释率（原代码缺失会报错）

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = group)) +
  geom_point(size = 3, alpha = 0.85) +
  # 手动指定分组颜色：Normal蓝色，Tumor红色
  scale_color_manual(values = c("Normal" = "#1f77b4", "Tumor" = "#d62728")) +
  labs(title = "PCA of TCGA-COAD samples",
       x = paste0("PC1 (", percent_var[1], "%)"),
       y = paste0("PC2 (", percent_var[2], "%)"), 
       color = "Group") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))

safe_ggsave("results/figures/DEG/PCA_plot.png", p_pca, 7, 6)

############################################################
## MA plot with ggplot2
############################################################

library(ggplot2)

deg_ma <- as.data.frame(res)
deg_ma$gene <- rownames(deg_ma)
deg_ma <- deg_ma[!is.na(deg_ma$padj) & !is.na(deg_ma$baseMean), ]

deg_ma$Group <- "Not significant"
deg_ma$Group[deg_ma$padj < 0.05 & deg_ma$log2FoldChange > 1] <- "Up"
deg_ma$Group[deg_ma$padj < 0.05 & deg_ma$log2FoldChange < -1] <- "Down"

p_ma <- ggplot(
  deg_ma,
  aes(x = log10(baseMean + 1), y = log2FoldChange, color = Group)
) +
  geom_point(alpha = 0.6, size = 0.8) +
  geom_hline(yintercept = c(-1, 1), linetype = 2) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_color_manual(
    values = c(
      "Up" = "#D73027",
      "Down" = "#4575B4",
      "Not significant" = "grey75"
    )
  ) +
  labs(
    title = "MA plot of TCGA-COAD differential expression",
    x = "log10 mean normalized expression",
    y = "log2 fold change"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank(),
    text = element_text(size = 13)
  )

ggsave(
  filename = "results/figures/DEG/MA_plot_ggplot.png",
  plot = p_ma,
  width = 7,
  height = 6,
  dpi = 600
)

ggsave(
  filename = "results/figures/DEG/MA_plot_ggplot.pdf",
  plot = p_ma,
  width = 7,
  height = 6
)

# Volcano
p_volcano <- ggplot(res_df, aes(log2FoldChange, -log10(padj), color = DEG_Group)) +
  geom_point(size = 1.1, alpha = 0.65) +
  geom_vline(xintercept = c(-1.5, 1.5), linetype = 2) +
  geom_hline(yintercept = -log10(0.01), linetype = 2) +
  scale_color_manual(values = c("Up" = "#D73027", "Down" = "#4575B4", "Not significant" = "grey70")) +
  labs(title = "Differentially expressed genes", x = "log2 fold change", y = "-log10(adj. P)", color = "") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))
safe_ggsave("results/figures/DEG/Volcano_plot.png", p_volcano, 8, 7)

# Top50 heatmap
symbol_labels <- res_df$gene_name
names(symbol_labels) <- res_df$GeneID
top50 <- res_df %>% arrange(padj) %>% slice_head(n = 50) %>% pull(GeneID)
heat_mat <- vsd_mat[top50, , drop = FALSE]
heat_mat <- t(scale(t(heat_mat)))
rownames(heat_mat) <- ifelse(is.na(symbol_labels[rownames(heat_mat)]), rownames(heat_mat), symbol_labels[rownames(heat_mat)])
ann_col <- data.frame(Group = coldata$group)
rownames(ann_col) <- rownames(coldata)
png("results/figures/DEG/Heatmap_Top50.png", width = 3000, height = 2400, res = 300)
pheatmap(heat_mat,
         annotation_col = ann_col,
         show_colnames = FALSE,
         show_rownames = TRUE,
         fontsize_row = 6,
         color = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
         border_color = NA,
         main = "Top 50 differentially expressed genes")
dev.off()

############################################################
## 3. WGCNA
############################################################

message("========== Step 3: WGCNA ==========")

# Use top variable genes for WGCNA to reduce noise and memory burden.
mad_values <- apply(vsd_mat, 1, mad)
max_wgcna_genes <- 8000
wgcna_genes <- names(sort(mad_values, decreasing = TRUE))[seq_len(min(max_wgcna_genes, length(mad_values)))]
wgcna_mat <- vsd_mat[wgcna_genes, , drop = FALSE]
datExpr0 <- t(wgcna_mat)

sampleTree <- hclust(dist(datExpr0), method = "average")
png("results/figures/WGCNA/Sample_Clustering.png", width = 2800, height = 1200, res = 300)
plot(sampleTree, main = "Sample clustering to detect outliers", xlab = "", sub = "", cex = 0.5)
dev.off()

# Conservative outlier handling: do not remove samples unless a user-specified cut is defensible.
# Here we keep all samples by default. If needed, inspect Sample_Clustering.png manually.
datExpr <- datExpr0
clinical_wgcna <- clinical[rownames(datExpr), , drop = FALSE]

gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  clinical_wgcna <- clinical_wgcna[rownames(datExpr), , drop = FALSE]
}

powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 5)
fit_df <- sft$fitIndices
candidate_power <- fit_df[fit_df$SFT.R.sq >= 0.85, ]
if (nrow(candidate_power) > 0) {
  softPower <- candidate_power$Power[1]
  selection_rule <- "minimum power with scale-free topology fit R^2 >= 0.85"
} else if (!is.na(sft$powerEstimate)) {
  softPower <- sft$powerEstimate
  selection_rule <- "WGCNA powerEstimate"
} else {
  softPower <- fit_df$Power[which.max(fit_df$SFT.R.sq)]
  selection_rule <- "power with maximum scale-free topology fit R^2"
}

save_csv(fit_df, "results/tables/WGCNA_SoftThreshold_FitIndices.csv")
save_csv(data.frame(selected_softPower = softPower, selection_rule = selection_rule),
         "results/tables/WGCNA_Selected_SoftPower.csv")

p_soft1 <- ggplot(fit_df, aes(Power, SFT.R.sq)) +
  geom_line() + geom_point() +
  geom_hline(yintercept = 0.85, linetype = 2, color = "red") +
  geom_vline(xintercept = softPower, linetype = 2) +
  theme_bw() + labs(title = "Scale-free topology fit", y = "Signed R^2")
p_soft2 <- ggplot(fit_df, aes(Power, mean.k.)) +
  geom_line() + geom_point() +
  geom_vline(xintercept = softPower, linetype = 2) +
  theme_bw() + labs(title = "Mean connectivity", y = "Mean connectivity")
# Save separately for robust portability.
safe_ggsave("results/figures/WGCNA/SoftThreshold_SFT.png", p_soft1, 6, 5)
safe_ggsave("results/figures/WGCNA/SoftThreshold_Connectivity.png", p_soft2, 6, 5)

net <- blockwiseModules(
  datExpr,
  power = softPower,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  saveTOMs = FALSE,
  verbose = 3
)

moduleColors <- labels2colors(net$colors)
MEs <- moduleEigengenes(datExpr, colors = moduleColors)$eigengenes
MEs <- orderMEs(MEs)

png("results/figures/WGCNA/Dendro_ModuleColors.png", width = 2600, height = 1400, res = 300)
plotDendroAndColors(net$dendrograms[[1]],
                    moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    main = "Gene dendrogram and module colors")
dev.off()

traitData <- data.frame(Tumor = ifelse(clinical_wgcna$group == "Tumor", 1, 0))
rownames(traitData) <- rownames(clinical_wgcna)
traitData <- traitData[match(rownames(datExpr), rownames(traitData)), , drop = FALSE]

moduleTraitCor <- cor(MEs, traitData, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))

save_csv(data.frame(Module = rownames(moduleTraitCor), Tumor_cor = moduleTraitCor[, 1], Tumor_pvalue = moduleTraitPvalue[, 1]),
         "results/tables/WGCNA_ModuleTrait.csv")

png("results/figures/WGCNA/Module_Trait_Relationship.png", width = 1800, height = max(1800, 90 * nrow(moduleTraitCor)), res = 300)
textMatrix <- paste0(sprintf("%.2f", moduleTraitCor), "\n(", format(moduleTraitPvalue, scientific = TRUE, digits = 2), ")")
dim(textMatrix) <- dim(moduleTraitCor)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(traitData),
               yLabels = rownames(moduleTraitCor),
               ySymbols = rownames(moduleTraitCor),
               colorLabels = FALSE,
               colors = blueWhiteRed(100),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1, 1),
               main = "Module-trait relationships")
dev.off()

keyME <- rownames(moduleTraitCor)[which.max(abs(moduleTraitCor[, "Tumor"]))]
keyModule <- gsub("^ME", "", keyME)
moduleGenes <- moduleColors == keyModule
moduleGeneIDs <- colnames(datExpr)[moduleGenes]

# GS and MM for key module
MM_all <- cor(datExpr, MEs, use = "p")
GS_Tumor <- cor(datExpr, traitData$Tumor, use = "p")
keyMM_col <- grep(paste0("^ME", keyModule, "$"), colnames(MM_all))
if (length(keyMM_col) != 1) stop("Cannot identify module eigengene column for key module: ", keyModule)

GSMM <- data.frame(
  GeneID = colnames(datExpr),
  ensembl_clean = clean_ensembl(colnames(datExpr)),
  Module = moduleColors,
  MM = MM_all[, keyMM_col],
  GS_Tumor = GS_Tumor[, 1]
)
GSMM <- dplyr::left_join(GSMM,
                         gene_anno %>% dplyr::select(ensembl_clean, gene_name, gene_type) %>% dplyr::distinct(),
                         by = "ensembl_clean")

CandidateHub <- GSMM %>%
  dplyr::filter(Module == keyModule, abs(MM) > 0.80, abs(GS_Tumor) > 0.30)

DEG_loose <- res_df %>% dplyr::filter(DEG_loose)
CoreHub <- CandidateHub %>%
  dplyr::inner_join(DEG_loose %>% dplyr::select(ensembl_clean, log2FoldChange, padj, DEG_Group), by = "ensembl_clean") %>%
  dplyr::arrange(padj, desc(abs(MM)), desc(abs(GS_Tumor)))

# If intersection is too small, keep CandidateHub as a fallback set for downstream scRNA localization.
HubForDownstream <- if (nrow(CoreHub) >= 3) CoreHub else CandidateHub
HubForDownstream <- HubForDownstream %>% dplyr::filter(!is.na(gene_name), gene_name != "")

save_csv(GSMM, paste0("results/tables/WGCNA_", keyModule, "_GS_MM_Table.csv"))
save_csv(CandidateHub, paste0("results/tables/WGCNA_", keyModule, "_CandidateHubGenes.csv"))
save_csv(CoreHub, paste0("results/tables/CoreHubGenes_DEG_and_WGCNA.csv"))
save_csv(HubForDownstream, paste0("results/tables/HubGenes_for_downstream.csv"))

p_gsmm <- ggplot(GSMM %>% dplyr::filter(Module == keyModule), aes(abs(MM), abs(GS_Tumor))) +
  geom_point(color = "grey65", alpha = 0.65, size = 1.4) +
  geom_point(data = CandidateHub, color = "red", size = 2.0) +
  geom_vline(xintercept = 0.80, linetype = 2) +
  geom_hline(yintercept = 0.30, linetype = 2) +
  theme_bw() +
  labs(title = paste0("GS-MM analysis of ", keyModule, " module"),
       x = "|Module membership|", y = "|Gene significance for Tumor|") +
  theme(plot.title = element_text(hjust = 0.5))
safe_ggsave(paste0("results/figures/WGCNA/", keyModule, "_GS_MM_Tumor.png"), p_gsmm, 7, 6)

save(datExpr, clinical_wgcna, moduleColors, MEs, moduleTraitCor, moduleTraitPvalue,
     GSMM, CandidateHub, CoreHub, HubForDownstream, keyModule,
     file = "results/RData/WGCNA_results.RData")

############################################################
## 4. GO / KEGG enrichment
############################################################

message("========== Step 4: GO / KEGG enrichment ==========")

enrich_symbols <- unique(na.omit(HubForDownstream$gene_name))
if (length(enrich_symbols) >= 3) {
  gene_df <- bitr(enrich_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  entrez_ids <- unique(gene_df$ENTREZID)
  if (length(entrez_ids) >= 3) {
    ego_bp <- enrichGO(entrez_ids, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
    ego_cc <- enrichGO(entrez_ids, OrgDb = org.Hs.eg.db, ont = "CC", readable = TRUE)
    ego_mf <- enrichGO(entrez_ids, OrgDb = org.Hs.eg.db, ont = "MF", readable = TRUE)
    ekegg <- enrichKEGG(entrez_ids, organism = "hsa")
    if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
      ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db)
    }
    save_csv(as.data.frame(ego_bp), "results/tables/GO_BP_enrichment.csv")
    save_csv(as.data.frame(ego_cc), "results/tables/GO_CC_enrichment.csv")
    save_csv(as.data.frame(ego_mf), "results/tables/GO_MF_enrichment.csv")
    save_csv(as.data.frame(ekegg), "results/tables/KEGG_enrichment.csv")
    if (nrow(as.data.frame(ego_bp)) > 0) {
      p_bp <- dotplot(ego_bp, showCategory = 15) + ggtitle("GO BP enrichment")
      safe_ggsave("results/figures/enrichment/GO_BP_dotplot.png", p_bp, 8, 6)
    }
    if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
      p_kegg <- dotplot(ekegg, showCategory = 15) + ggtitle("KEGG enrichment")
      safe_ggsave("results/figures/enrichment/KEGG_dotplot.png", p_kegg, 8, 6)
    }
  }
} else {
  writeLines("Fewer than 3 hub gene symbols are available for enrichment.",
             "results/tables/Enrichment_Not_Run.txt")
}

############################################################
## 5. Survival and diagnostic ROC analyses
############################################################

message("========== Step 5: Survival and ROC ==========")

# Survival should be performed on tumor samples only and patient-level expression.
# Diagnostic ROC compares Tumor vs Normal expression, not survival status.

hub_symbols <- unique(na.omit(HubForDownstream$gene_name))
expr_symbol_vst <- collapse_by_symbol(vsd_mat, gene_anno, stat = "mean")

# Diagnostic ROC: Tumor vs Normal

dir.create("results/figures/survival", recursive = TRUE, showWarnings = FALSE)

roc_results <- list()

for (g in hub_symbols[hub_symbols %in% rownames(expr_symbol_vst)]) {
  
  df <- data.frame(
    expr = as.numeric(expr_symbol_vst[g, ]),
    group = coldata$group
  )
  
  df <- df[!is.na(df$expr) & !is.na(df$group), ]
  df$group <- factor(df$group, levels = c("Normal", "Tumor"))
  
  if (length(unique(df$group)) == 2) {
    
    roc_obj <- pROC::roc(
      response = df$group,
      predictor = df$expr,
      levels = c("Normal", "Tumor"),
      quiet = TRUE
    )
    
    roc_results[[g]] <- data.frame(
      gene_name = g,
      AUC = as.numeric(pROC::auc(roc_obj))
    )
    
    png(
      filename = paste0("results/figures/survival/ROC_", g, ".png"),
      width = 1600,
      height = 1400,
      res = 300
    )
    
    pROC::plot.roc(
      roc_obj,
      col = "red",
      legacy.axes = TRUE,
      main = paste0("Diagnostic ROC: ", g)
    )
    
    dev.off()
  }
}
if (length(roc_results) > 0) save_csv(dplyr::bind_rows(roc_results), "results/tables/Diagnostic_ROC_AUC.csv")

# Survival: use TCGA clinical fields where available.
clin_surv <- clinical[clinical$group == "Tumor", , drop = FALSE]
clin_surv$days_to_death_num <- suppressWarnings(as.numeric(clin_surv$days_to_death))
clin_surv$days_to_last_follow_up_num <- suppressWarnings(as.numeric(clin_surv$days_to_last_follow_up))
clin_surv$time <- ifelse(!is.na(clin_surv$days_to_death_num),
                         clin_surv$days_to_death_num,
                         clin_surv$days_to_last_follow_up_num)
clin_surv$status <- ifelse(clin_surv$vital_status == "Dead", 1, 0)
clin_surv <- clin_surv[!is.na(clin_surv$time) & clin_surv$time > 0 & !is.na(clin_surv$status), , drop = FALSE]

tumor_samples <- rownames(clin_surv)
expr_tumor <- expr_symbol_vst[, intersect(colnames(expr_symbol_vst), tumor_samples), drop = FALSE]
clin_surv <- clin_surv[colnames(expr_tumor), , drop = FALSE]

surv_summary <- list()
if (nrow(clin_surv) >= 30) {
  for (g in hub_symbols[hub_symbols %in% rownames(expr_tumor)]) {
    surv_df <- data.frame(
      time = clin_surv$time,
      status = clin_surv$status,
      expr = as.numeric(expr_tumor[g, ]),
      sample_barcode = rownames(clin_surv)
    )
    surv_df <- surv_df[complete.cases(surv_df), ]
    if (nrow(surv_df) < 30 || length(unique(surv_df$status)) < 2) next
    surv_df$expr_group <- ifelse(surv_df$expr > median(surv_df$expr, na.rm = TRUE), "High", "Low")
    fit <- survfit(Surv(time, status) ~ expr_group, data = surv_df)
    p <- ggsurvplot(fit, data = surv_df, pval = TRUE, risk.table = TRUE,
                    title = paste0("Overall survival by ", g, " expression"))
    ggsave(paste0("results/figures/survival/KM_", g, ".png"), plot = p$plot, width = 7, height = 6, dpi = 600)
    sdif <- survdiff(Surv(time, status) ~ expr_group, data = surv_df)
    pval <- 1 - pchisq(sdif$chisq, length(sdif$n) - 1)
    surv_summary[[g]] <- data.frame(gene_name = g, logrank_p = pval, n = nrow(surv_df))
  }
}
if (length(surv_summary) > 0) save_csv(dplyr::bind_rows(surv_summary), "results/tables/Survival_Logrank_Summary.csv")

############################################################
## 6. scRNA-seq analysis and Hub gene cell localization
############################################################

message("========== Step 6: scRNA-seq and cell localization ==========")

# GSE132465 is large. This function tries to create a Seurat object from common GEO supplementary formats.
# If the downloaded supplementary files are not a directly readable matrix, place a processed matrix or Seurat RDS
# under data/scRNA/ and re-run this section.

read_first_expression_matrix <- function(sc_dir = "data/scRNA", geo_id = "GSE132465") {
  seurat_rds <- list.files(sc_dir, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
  if (length(seurat_rds) > 0) {
    obj <- readRDS(seurat_rds[1])
    if (inherits(obj, "Seurat")) return(obj)
  }

  if (length(list.files(sc_dir, recursive = TRUE)) == 0) {
    GEOquery::getGEOSuppFiles(geo_id, baseDir = sc_dir, makeDirectory = TRUE)
  }

  files <- list.files(sc_dir, recursive = TRUE, full.names = TRUE)

  # 10x h5
  h5_files <- files[grepl("\\.h5$", files, ignore.case = TRUE)]
  if (length(h5_files) > 0) {
    mat <- Read10X_h5(h5_files[1])
    if (is.list(mat)) mat <- mat[[1]]
    return(CreateSeuratObject(counts = mat, project = geo_id))
  }

  # 10x folder with matrix.mtx / barcodes / features or genes
  mtx_files <- files[grepl("matrix\\.mtx(\\.gz)?$", files, ignore.case = TRUE)]
  if (length(mtx_files) > 0) {
    tenx_dir <- dirname(mtx_files[1])
    mat <- Read10X(data.dir = tenx_dir)
    if (is.list(mat)) mat <- mat[[1]]
    return(CreateSeuratObject(counts = mat, project = geo_id))
  }

  # Wide expression matrix. Use the first likely file.
  tab_files <- files[grepl("\\.(txt|tsv|csv)(\\.gz)?$", files, ignore.case = TRUE)]
  if (length(tab_files) > 0) {
    for (f in tab_files) {
      message("Trying to read scRNA matrix: ", f)
      dt <- tryCatch(data.table::fread(f), error = function(e) NULL)
      if (is.null(dt) || ncol(dt) < 10 || nrow(dt) < 100) next
      gene_col <- 1
      gene_names <- dt[[gene_col]]
      expr_dt <- dt[, -gene_col, with = FALSE]
      expr_mat <- as.matrix(expr_dt)
      mode(expr_mat) <- "numeric"
      rownames(expr_mat) <- make.unique(as.character(gene_names))
      colnames(expr_mat) <- make.unique(colnames(expr_dt))
      return(CreateSeuratObject(counts = Matrix(expr_mat, sparse = TRUE), project = geo_id))
    }
  }

  stop("No readable scRNA expression matrix was found under ", sc_dir,
       ". Please download a processed count matrix or Seurat object manually.")
}

sc <- tryCatch(read_first_expression_matrix("data/scRNA", "GSE132465"), error = function(e) {
  writeLines(conditionMessage(e), "results/tables/scRNA_Loading_Failed.txt")
  NULL
})

if (!is.null(sc)) {
  sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
  sc <- subset(sc, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
  sc <- NormalizeData(sc)
  sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 3000)
  sc <- ScaleData(sc, verbose = FALSE)
  sc <- RunPCA(sc, npcs = 30, verbose = FALSE)
  sc <- FindNeighbors(sc, dims = 1:20)
  sc <- FindClusters(sc, resolution = 0.5)
  sc <- RunUMAP(sc, dims = 1:20)

  p_umap_cluster <- DimPlot(sc, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
    ggtitle("scRNA-seq clusters")
  safe_ggsave("results/figures/scRNA/UMAP_clusters.png", p_umap_cluster, 7, 6)

  # Marker-based coarse cell-type localization. If cell-type metadata exists, use it; otherwise assign by marker scores.
  meta_cols <- colnames(sc@meta.data)
  celltype_col <- meta_cols[grepl("cell.?type|annotation|celltype", meta_cols, ignore.case = TRUE)][1]

  if (!is.na(celltype_col)) {
    sc$cell_type_final <- sc@meta.data[[celltype_col]]
  } else {
    marker_list <- list(
      Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT20"),
      T_NK = c("CD3D", "CD3E", "TRAC", "NKG7", "GNLY"),
      B_cell = c("MS4A1", "CD79A", "CD79B"),
      Plasma = c("MZB1", "JCHAIN", "XBP1"),
      Myeloid = c("LYZ", "LST1", "S100A8", "S100A9"),
      Fibroblast = c("COL1A1", "COL1A2", "DCN", "LUM"),
      Endothelial = c("PECAM1", "VWF", "KDR"),
      Mast = c("TPSAB1", "TPSB2", "CPA3")
    )
    marker_list <- lapply(marker_list, function(x) intersect(x, rownames(sc)))
    marker_list <- marker_list[lengths(marker_list) > 0]
    sc <- AddModuleScore(sc, features = marker_list, name = "CTScore")
    score_cols <- paste0("CTScore", seq_along(marker_list))
    score_mat <- sc@meta.data[, score_cols, drop = FALSE]
    sc$cell_type_final <- names(marker_list)[max.col(score_mat, ties.method = "first")]

    cluster_celltype <- sc@meta.data %>%
      dplyr::count(seurat_clusters, cell_type_final) %>%
      dplyr::group_by(seurat_clusters) %>%
      dplyr::slice_max(n, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()
    save_csv(cluster_celltype, "results/tables/scRNA_Cluster_CellType_Assignment.csv")
  }

  p_umap_celltype <- DimPlot(sc, reduction = "umap", group.by = "cell_type_final", label = TRUE, repel = TRUE) +
    ggtitle("scRNA-seq coarse cell types")
  safe_ggsave("results/figures/scRNA/UMAP_cell_types.png", p_umap_celltype, 8, 6)

  sc_hub_genes <- intersect(hub_symbols, rownames(sc))
  save_csv(data.frame(hub_gene = hub_symbols, found_in_scRNA = hub_symbols %in% rownames(sc)),
           "results/tables/scRNA_HubGene_Found_Missing.csv")

  if (length(sc_hub_genes) > 0) {
    # Feature plots for up to 12 genes.
    fp_genes <- head(sc_hub_genes, 12)
    p_feat <- FeaturePlot(sc, features = fp_genes, reduction = "umap", ncol = 3)
    safe_ggsave("results/figures/scRNA/HubGene_FeaturePlot.png", p_feat, 12, 10)

    p_dot <- DotPlot(sc, features = sc_hub_genes, group.by = "cell_type_final") +
      RotatedAxis() + ggtitle("Hub gene expression across cell types")
    safe_ggsave("results/figures/scRNA/HubGene_CellType_DotPlot.png", p_dot, max(8, length(sc_hub_genes) * 0.55), 5)

    p_vln <- VlnPlot(sc, features = fp_genes, group.by = "cell_type_final", pt.size = 0, ncol = 3)
    safe_ggsave("results/figures/scRNA/HubGene_CellType_Violin.png", p_vln, 12, 10)

    sc <- AddModuleScore(sc, features = list(sc_hub_genes), name = "HubScore")
    p_score_feat <- FeaturePlot(sc, features = "HubScore1", reduction = "umap") + ggtitle("Hub gene set score")
    safe_ggsave("results/figures/scRNA/HubScore_FeaturePlot.png", p_score_feat, 7, 6)
    p_score_vln <- VlnPlot(sc, features = "HubScore1", group.by = "cell_type_final", pt.size = 0) +
      ggtitle("Hub gene set score by cell type")
    safe_ggsave("results/figures/scRNA/HubScore_CellType_Violin.png", p_score_vln, 9, 5)

    avg_expr <- AverageExpression(sc, features = sc_hub_genes, group.by = "cell_type_final", assays = "RNA", slot = "data")$RNA
    avg_df <- as.data.frame(as.matrix(avg_expr)) %>%
      tibble::rownames_to_column("gene_name")
    save_csv(avg_df, "results/tables/scRNA_HubGene_CellType_AverageExpression.csv")

    det_list <- lapply(sc_hub_genes, function(g) {
      expr_vec <- FetchData(sc, vars = g)[, 1]
      data.frame(gene_name = g, cell_type = sc$cell_type_final, detected = expr_vec > 0) %>%
        dplyr::group_by(gene_name, cell_type) %>%
        dplyr::summarise(detection_fraction = mean(detected), .groups = "drop")
    })
    det_df <- dplyr::bind_rows(det_list)
    save_csv(det_df, "results/tables/scRNA_HubGene_CellType_DetectionFraction.csv")

    loc_df <- det_df %>%
      dplyr::left_join(
        avg_df %>% tidyr::pivot_longer(-gene_name, names_to = "cell_type", values_to = "average_expression"),
        by = c("gene_name", "cell_type")
      ) %>%
      dplyr::group_by(gene_name) %>%
      dplyr::arrange(desc(average_expression), desc(detection_fraction), .by_group = TRUE) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::rename(top_localized_cell_type = cell_type)
    save_csv(loc_df, "results/tables/scRNA_HubGene_CellType_Localization.csv")
  } else {
    writeLines("No bulk hub gene symbols were found in the scRNA-seq object rownames.",
               "results/tables/scRNA_HubGene_Localization_Not_Run.txt")
  }

  saveRDS(sc, "data/scRNA/GSE132465_Seurat_processed.rds")
}

############################################################
## 7. Final summary
############################################################

message("========== Step 7: Final summary ==========")
final_summary <- data.frame(
  TP_samples = sum(clinical$group == "Tumor"),
  NT_samples = sum(clinical$group == "Normal"),
  genes_before_filter = before_filter_genes,
  genes_after_filter = after_filter_genes,
  DEG_strict_up = DEG_summary$strict_up,
  DEG_strict_down = DEG_summary$strict_down,
  DEG_loose_up = DEG_summary$loose_up,
  DEG_loose_down = DEG_summary$loose_down,
  WGCNA_input_genes = ncol(datExpr),
  selected_softPower = softPower,
  keyModule = keyModule,
  keyModule_Tumor_cor = moduleTraitCor[keyME, "Tumor"],
  keyModule_Tumor_pvalue = moduleTraitPvalue[keyME, "Tumor"],
  CandidateHubGene_count = nrow(CandidateHub),
  CoreHubGene_DEG_WGCNA_count = nrow(CoreHub),
  DownstreamHubGene_count = nrow(HubForDownstream),
  scRNA_loaded = !is.null(sc)
)
save_csv(final_summary, "results/tables/Final_Result_Summary.csv")

writeLines(capture.output(sessionInfo()), "results/tables/sessionInfo.txt")
message("Pipeline finished. Key outputs are under results/tables and results/figures.")
