# Author: Inho Park
# Purpose: Run bulkSignalR for a visium dataset

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
plot_dir <- file.path(output_dir, "bulkSignalR", "figures")

n.proc <- 8 # number of cores to use
cl <- makeCluster(n.proc)
registerDoParallel(cl)
#stopCluster(cl)

args <- commandArgs(trailingOnly = TRUE)
sampleId <- args[1]
message("Run Sample: ", sampleId)

# sampleId
#sampleId <- "YS15-33091"
seurat_processed <- file.path(proxy_dir, paste0(sampleId, "_seurat_processed.RDS"))
seurat_obj <- readRDS(seurat_processed)

gene_expr_df <-  seurat_obj@assays$Spatial@counts
gene_expr_df <- as.data.frame(gene_expr_df)

bsrdm <- prepareDataset(counts=gene_expr_df, min.count = 1)
bsrdm <- learnParameters(bsrdm, min.positive = 2, plot.folder = plot_dir, filename = sampleId)

bsrinf <- initialInference(bsrdm, min.cor = -1)
stopCluster(cl)

saveRDS(bsrdm, file.path(proxy_dir, paste0(sampleId, "_bsrdm.RDS")))
saveRDS(bsrinf, file.path(proxy_dir, paste0(sampleId, "_bsrinf.RDS")))

#LRinter.dataframe <- LRinter(bsrinf)
#head(LRinter.dataframe)
#LRinter.dataframe.selected <- LRinter.dataframe[LRinter.dataframe$qval<= 0.001,]
#LRinter.dataframe.selected <- LRinter.dataframe.selected[order(LRinter.dataframe.selected$qval),]
#bsrinf.redBP <- reduceToBestPathway(bsrinf)
#bsrinf.L    <- reduceToLigand(bsrinf)
#bsrinf.R    <- reduceToReceptor(bsrinf)
#bsrinf.redP  <- reduceToPathway(bsrinf)  
#bsrinf.redPBP <- reduceToBestPathway(bsrinf.redP) 
#bsrsig.redBP <- getLRGeneSignatures(bsrinf.redBP, qval.thres=0.001)
#scoresLR <- scoreLRGeneSignatures(bsrdm, bsrsig.redBP, name.by.pathway=FALSE)
#simpleHeatmap(scoresLR[1:20,], path=plot_dir, filename=paste0("/",sampleId, "_scoresLR"), column.names=TRUE, height=5, width=9, pointsize=10, hcl.palette="Cividis")
#bsrsig.redPBP <- getLRGeneSignatures(bsrinf.redPBP, qval.thres=0.001)
#     
#scoresPathway <- scoreLRGeneSignatures(bsrdm, bsrsig.redPBP, name.by.pathway=TRUE)    
#simpleHeatmap(scoresPathway[1:10,], 
#                   path = plot_dir,
#                   filename = paste0("/",sampleId, "_scoresPathway"),
#                   column.names = TRUE, 
#                   width = 9, 
#                   height = 5, 
#                   pointsize = 12,
#                   hcl.palette = "Blue-Red 2"
#                   )
#pdf(file.path(plot_dir, paste0(sampleId, "_scoresPathway.pdf")), width=20, height=8)
#pheatmap(scoresPathway[1:10,], fontsize_col=7)
#dev.off()
