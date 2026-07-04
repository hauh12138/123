rm(list = ls())
gc()
options(stringsAsFactors = FALSE)
options(timeout = 600)
setwd("E:\\HuaweiMoveData\\Users\\潘懿\\Desktop\\abab")
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)

library(WGCNA)

library(clusterProfiler)
library(org.Hs.eg.db)

library(ggplot2)
library(pheatmap)

library(survival)
library(survminer)
library(pROC)

library(Seurat)
library(GEOquery)
library(dplyr)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

dir.create("data",
           showWarnings = FALSE)
dir.create("data/raw",
           showWarnings = FALSE)
dir.create("data/processed",
           showWarnings = FALSE)
dir.create("data/clinical",
           showWarnings = FALSE)
dir.create("results",
           showWarnings = FALSE)
dir.create("results/figures",
           showWarnings = FALSE)
dir.create("results/tables",
           showWarnings = FALSE)

query <- GDCquery(
  project="TCGA-COAD",
  data.category="Transcriptome Profiling",
  data.type="Gene Expression Quantification",
  workflow.type="STAR - Counts"
)

GDCdownload(query)
data <- GDCprepare(query)

expr <- assay(data)
clinical <- colData(data)

group <- ifelse(substr(colnames(expr),14,15)=="11","Normal","Tumor")
group <- factor(group, levels=c("Normal","Tumor"))

expr <- expr[rowSums(expr > 10) > 5, ]

coldata <- data.frame(row.names=colnames(expr), group=group)

dds <- DESeqDataSetFromMatrix(expr, coldata, design=~group)
dds <- dds[rowSums(counts(dds))>10,]

dds <- DESeq(dds)
res <- results(dds)
res <- na.omit(res)

deg <- as.data.frame(res)
deg$gene <- rownames(deg)

vsd <- vst(dds)
mat <- assay(vsd)

pca <- prcomp(t(mat))

ggplot(data.frame(pca$x), aes(PC1,PC2,color=group))+
  geom_point(size=3)+
  theme_classic()

# =========================
# WGCNA + sample trimming（SCI标准版）
# =========================

datExpr0 <- t(mat)

# 1. sample clustering
sampleTree <- hclust(dist(datExpr0), method="average")

plot(sampleTree,
     main="Sample clustering to detect outliers",
     xlab="", sub="")

# 2. cut tree (outlier detection)
cutHeight <- 100

clust <- cutreeStaticTree(sampleTree,
                          cutHeight = cutHeight,
                          minSize = 10)

keepSamples <- (clust == 1)

datExpr0 <- datExpr0[keepSamples, ]

# 3. gene filtering
gsg <- goodSamplesGenes(datExpr0)
datExpr <- datExpr0[gsg$goodSamples, gsg$goodGenes]

# 4. soft threshold
powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector=powers)
softPower <- sft$powerEstimate

# 5. network
net <- blockwiseModules(
  datExpr,
  power=softPower,
  TOMType="unsigned",
  minModuleSize=30,
  mergeCutHeight=0.25,
  numericLabels=TRUE,
  saveTOMs=TRUE,
  saveTOMFileBase="TOM"
)

moduleColors <- labels2colors(net$colors)
MEs <- net$MEs

plotDendroAndColors(net$dendrograms[[1]],
                    moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    main="Gene dendrogram")

trait <- data.frame(Tumor=as.numeric(group=="Tumor"))
rownames(trait) <- rownames(datExpr)

moduleTraitCor <- cor(MEs, trait, use="p")

pheatmap(moduleTraitCor,
         color=colorRampPalette(c("blue","white","red"))(50),
         main="Module-trait relationships")

geneModuleMembership <- cor(datExpr, MEs, use="p")
geneTraitSignificance <- cor(datExpr, trait$Tumor)

hub_genes <- names(which(
  abs(geneModuleMembership[,1]) > 0.8 &
    abs(geneTraitSignificance) > 0.2
))

gene_df <- bitr(hub_genes,
                fromType="SYMBOL",
                toType="ENTREZID",
                OrgDb=org.Hs.eg.db)

genes <- gene_df$ENTREZID

ego_bp <- enrichGO(genes,
                   OrgDb=org.Hs.eg.db,
                   ont="BP",
                   readable=TRUE)

ego_cc <- enrichGO(genes,
                   OrgDb=org.Hs.eg.db,
                   ont="CC")

ego_mf <- enrichGO(genes,
                   OrgDb=org.Hs.eg.db,
                   ont="MF")

ekegg <- enrichKEGG(genes,
                    organism="hsa")

ekegg <- setReadable(ekegg, OrgDb=org.Hs.eg.db)

dotplot(ego_bp, showCategory=15)
dotplot(ekegg, showCategory=15)
cnetplot(ego_bp, showCategory=5)

clinical <- as.data.frame(clinical)

clinical$time <- ifelse(
  is.na(clinical$days_to_death),
  clinical$days_to_last_follow_up,
  clinical$days_to_death
)

clinical$status <- ifelse(clinical$vital_status=="Dead",1,0)

colnames(expr) <- substr(colnames(expr),1,16)
rownames(clinical) <- substr(rownames(clinical),1,16)

gene <- hub_genes[1]

gene_expr <- as.data.frame(t(expr[gene,]))
colnames(gene_expr) <- "expr"

surv_data <- cbind(clinical, gene_expr)
surv_data <- na.omit(surv_data)

surv_data$group <- ifelse(surv_data$expr > median(surv_data$expr),"High","Low")

fit <- survfit(Surv(time,status)~group,data=surv_data)

ggsurvplot(fit,data=surv_data,pval=TRUE)

roc_obj <- roc(surv_data$status, surv_data$expr)
plot(roc_obj, col="red")
auc(roc_obj)

geo <- getGEO("GSE132465",GSEMatrix=TRUE)
sce <- geo[[1]]

sc <- CreateSeuratObject(exprs(sce))

# QC
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern="^MT-")

sc <- subset(sc,
             subset = nFeature_RNA > 200 &
               nFeature_RNA < 5000 &
               percent.mt < 10)

# standard pipeline
sc <- NormalizeData(sc)
sc <- FindVariableFeatures(sc)
sc <- ScaleData(sc)

sc <- RunPCA(sc)
ElbowPlot(sc)

sc <- FindNeighbors(sc, dims=1:20)
sc <- FindClusters(sc, resolution=0.5)

sc <- RunUMAP(sc, dims=1:20)

DimPlot(sc, label=TRUE)

# hub gene validation
FeaturePlot(sc, features=hub_genes[1:3])

# module score
sc <- AddModuleScore(sc,
                     features=list(hub_genes),
                     name="HubScore")

FeaturePlot(sc, features="HubScore1")
