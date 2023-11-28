library(config)
library(BulkSignalR)
library(igraph)
library(Seurat)
library(dplyr)
library(doParallel)
library(pheatmap)
library(ggpubr)
library(ggplot2)

# Load config data
config <- config::get(file="config.yaml")
proxy_dir <- config$proxy
output_dir <- config$output
data_dir <- config$data
plot_dir <- file.path(output_dir, "bulkSignalR", "figures")


args <- commandArgs(trailingOnly = TRUE)
sampleId <- args[1]
message("Run Sample: ", sampleId)

#sampleId = "YS15-33091"
bsrdm <- readRDS(file.path(proxy_dir, paste0(sampleId, "_bsrdm.RDS")))
bsrinf <- readRDS(file.path(proxy_dir, paste0(sampleId, "_bsrinf.RDS")))

bsrinf.L    <- reduceToLigand(bsrinf)
bsrinf.R    <- reduceToReceptor(bsrinf)
bsrinf.redBP <- reduceToBestPathway(bsrinf)

bsrinf.redP  <- reduceToPathway(bsrinf)  
bsrinf.redPBP <- reduceToBestPathway(bsrinf.redP) 

bsrsig.redBP <- getLRGeneSignatures(bsrinf.redBP, qval.thres=0.001)
scoresLR <- scoreLRGeneSignatures(bsrdm, bsrsig.redBP, name.by.pathway=FALSE)

bsrsig.redPBP <- getLRGeneSignatures(bsrinf.redPBP, qval.thres=0.001)    
scoresPathway <- scoreLRGeneSignatures(bsrdm, bsrsig.redPBP, name.by.pathway=TRUE)

names(scoresPathway[1,][order(-scoresPathway[1,])][1:100]
top100 <- names(scoresPathway[1,][order(-scoresPathway[1,])][1:100])

# reduce to the clusters
seurat_processed <- file.path(proxy_dir, paste0(sampleId, "_seurat_processed.RDS"))
seurat_obj <- readRDS(seurat_processed)

patho_df <- data.frame(seurat_obj$patho)
patho_df$barcode <- rownames(patho_df)
colnames(patho_df) <- c("PATHO", "barcode")

patho_df <- patho_df[complete.cases(patho_df),]
patho_categories <- unique(patho_df$PATHO)

all_scores <- data.frame()
for (patho_category in patho_categories) {
    cell_selector <- patho_df$PATHO==patho_category    
    cell_ids = patho_df[cell_selector,]$barcode

    scoresPathway_categories = scoresPathway[,cell_ids]
    group_mean <- apply(scoresPathway_categories, 1, mean)
    group_mean <- t(data.frame(group_mean))
    #colnames(group_mean) <- patho_category

    #all_scores <- cbind(all_scores, group_mean)
    #print(group_mean)
    all_scores <- rbind(all_scores, group_mean)    
    #kk <- apply(scoresPathway_categories, 2, mean)
    #print(kk)
    #spots_groups[[patho]] <- cell_ids
}
rownames(all_scores) <- patho_categories
all_scores <- t(all_scores)

pdf(file.path(plot_dir, paste0(sampleId, "_patho_pathways.pdf")), width=7.5, height=10)
pheatmap(all_scores)
dev.off()
