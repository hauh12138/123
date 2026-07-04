# 123
12345
# TCGA-COAD scRNA Hub Gene Analysis

## Project title

基于 TCGA 与单细胞转录组数据筛选结肠腺癌核心 Hub 基因及细胞定位研究

## Main scripts

- `scripts/111.R`: original analysis script
- `scripts/COAD_TCGA_scRNA_Hub_pipeline_revised.R`: revised analysis pipeline

## Main result folders

- `results/figures/`: generated plots
- `results/tables/`: result tables

## Analysis workflow

1. TCGA-COAD STAR-Counts download
2. DESeq2 differential expression analysis
3. WGCNA module analysis
4. Hub gene selection
5. GO/KEGG enrichment
6. Survival and ROC analysis
7. scRNA-seq cell localization analysis

## Notes

Large raw data files are not included. The analysis can be reproduced by running the scripts in order.
