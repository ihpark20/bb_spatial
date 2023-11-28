# Author: Inho Park
# Purpose: Run bulkSignalR for a visium dataset

library(config)
library(BulkSignalR)
library(igraph)
library(Seurat)
library(dplyr)

# Load config data
config <- config::get(file="config.yaml")
proxy_dir <- config$proxy
output_dir <- config$output
plot_dir <- file.path(output_dir, "bulkSignalR", "figures")

#plot_folder <- file.path(pr)
# sampleId
sampleId <- "YS15-33091"
seurat_processed <- file.path(proxy_dir, paste0(sampleId, "_seurat_processed.RDS"))
seurat_obj <- readRDS(seurat_processed)

# create groups of spots
# 최초 clustering의 resolution을 좀 더 높게 설정하는 것도 생각해봐야할 듯
all_clusters = unique(seurat_obj@meta.data$seurat_clusters)
spots_groups = list()

for (cluster in all_clusters) {
    cell_selector <- seurat_obj$seurat_clusters==cluster    
    cell_ids = colnames(seurat_obj)[cell_selector]
    spots_groups[[cluster]] <- cell_ids
    print(length(cell_ids))
}

# show the first 10 spots in each cluster
spots_groups$`1`[1:10]
spots_groups$`2`[1:10]

# aggregate expression of spots in each cluster (raw counts)
aggregateCount <- function(seurat_data, assay_name, spot_names_in_cluster) {
    aggregated_gene_count <- apply(GetAssay(seurat_data, assay=assay_name)[, spot_names_in_cluster], 1, sum)
    aggregated_gene_count
}

#g1 <- aggregateCount(seurat_obj, "Spatial", spots_groups$`1`)
aggregated_list <- lapply(spots_groups, function(x) aggregateCount(seurat_obj, "Spatial", x))
gene_expr_df <- as.data.frame(aggregated_list)

bsrdm <- prepareDataset(counts=gene_expr_df)
bsrdm <- learnParameters(bsrdm, plot.folder = plot_dir, filename = sampleId)

bsrinf <- initialInference(bsrdm)
LRinter.dataframe <- LRinter(bsrinf)

head(LRinter.dataframe)

LRinter.dataframe.selected <- LRinter.dataframe[LRinter.dataframe$qval<= 0.001,]
LRinter.dataframe.selected <- LRinter.dataframe.selected[order(LRinter.dataframe.selected$qval),]

bsrinf.redBP <- reduceToBestPathway(bsrinf)
bsrinf.L    <- reduceToLigand(bsrinf)
bsrinf.R    <- reduceToReceptor(bsrinf)

bsrinf.redP  <- reduceToPathway(bsrinf)  
bsrinf.redPBP <- reduceToBestPathway(bsrinf.redP) 

bsrsig.redBP <- getLRGeneSignatures(bsrinf.redBP, qval.thres=0.001)

scoresLR <- scoreLRGeneSignatures(bsrdm, bsrsig.redBP, name.by.pathway=FALSE)
simpleHeatmap(scoresLR[1:20,], path=plot_dir, filename=paste0("/",sampleId, "_scoresLR"), column.names=TRUE, height=5, width=9, pointsize=10, hcl.palette="Cividis")

bsrsig.redPBP <- getLRGeneSignatures(bsrinf.redPBP, qval.thres=0.001)
     
scoresPathway <- scoreLRGeneSignatures(bsrdm, bsrsig.redPBP, name.by.pathway=TRUE)
     
simpleHeatmap(scoresPathway[1:10,], 
                   path = plot_dir,
                   filename = paste0("/",sampleId, "_scoresPathway"),
                   column.names = TRUE, 
                   width = 9, 
                   height = 5, 
                   pointsize = 12,
                   hcl.palette = "Blue-Red 2"
                   )
