# Author: Inho Park
# Purpose: Run bulkSignalR for a visium dataset

library(config)
library(BulkSignalR)
library(igraph)
library(Seurat)
library(dplyr)
library(doParallel)

# Load config data
config <- config::get(file="config.yaml")
proxy_dir <- config$proxy
output_dir <- config$output
plot_dir <- file.path(output_dir, "bulkSignalR", "figures")


n.proc <- 8 # number of cores to use
cl <- makeCluster(n.proc)
registerDoParallel(cl)

# sampleId
sampleId <- "YS15-33091"
seurat_processed <- file.path(proxy_dir, paste0(sampleId, "_seurat_processed.RDS"))
seurat_obj <- readRDS(seurat_processed)

cells <- c("B.cells", "CAFs", "Cancer.Epithelial", "Endothelial", "Myeloid", "Normal.Epithelial", "PVL", "Plasmablasts", "T.cells")
cells_ratio <- paste0(cells, "_ratio")

major_cells_idx <- apply(seurat_obj@meta.data[,cells_ratio], 1, function(x) which(x==max(x)))
major_cells <- cells[major_cells_idx]
seurat_obj$major_cell <- major_cells



tissue_positions_list_csv_file <- file.path(data_dir, "01.raw", sampleId, "Spatial", "tissue_positions_list.csv")
scalefactors_json <- file.path(data_dir, "01.raw", sampleId, "Spatial", "scalefactors_json.json")

# from utils.R
tissue_positions_df <- get_tissue_positions(tissue_positions_list_csv_file)
spot_dist_matrix <- get_spot_distance_matrix(tissue_positions_df, scalefactors_json)

# define your interest gruop
tumor_annotations <- c("DCIS", "invasive carcinoma", "tumor", "cancer")
tumor_spots <- names(seurat_obj$patho[seurat_obj$patho %in% tumor_annotations])
non_tumor_spots <- names(seurat_obj$patho[!seurat_obj$patho %in% tumor_annotations])

# define tumor core spots, tumor periphery spots
distances_from_non_tumor <- get_min_distance_from_source_spots_to_target_spots(spot_dist_matrix, tumor_spots, non_tumor_spots)
tumor_core_spots <- names(distances_from_non_tumor[distances_from_non_tumor > 120])
tumor_periphery_spots <- names(distances_from_non_tumor[distances_from_non_tumor <= 120])

# define non-tumor boundary spots
distances_from_tumor <- get_min_distance_from_source_spots_to_target_spots(spot_dist_matrix, non_tumor_spots, tumor_spots)
tumor_near_spots <- names(distances_from_tumor[distances_from_tumor <= 120])
tumor_dist_spots <- names(distances_from_tumor[distances_from_tumor > 120])

seurat_obj$SpotType <- NA
seurat_obj@meta.data[tumor_near_spots,]$SpotType <- "TumorNear"
seurat_obj@meta.data[tumor_dist_spots,]$SpotType <- "TumorDist"
seurat_obj@meta.data[tumor_core_spots,]$SpotType <- "TumorCore"
seurat_obj@meta.data[tumor_periphery_spots,]$SpotType <- "TumorPeri"





# create groups of spots
# 최초 clustering의 resolution을 좀 더 높게 설정하는 것도 생각해봐야할 듯
all_clusters = unique(seurat_obj@meta.data$seurat_clusters)
all_pathos = unique(seurat_obj$patho)

min_spots = 5
spots_groups = list()
for (cluster in all_clusters) {
    for(patho in all_pathos) {
        cell_selector <- seurat_obj$seurat_clusters==cluster & seurat_obj$patho==patho
        if (sum(cell_selector, na.rm=TRUE) < min_spots) {next}
        #message(cluster, ":", patho)
        #print(sum(cell_selector, na.rm=TRUE))
        cell_ids = colnames(seurat_obj)[cell_selector]
        cell_ids = na.omit(cell_ids)

        spots_groups_names <- str_replace(paste0(cluster, "_", patho), "-", "_")
        spots_groups_names <- str_replace(spots_groups_names, " ", "_")

        spots_groups[[spots_groups_names]] <- cell_ids
    }
}

# show the first 10 spots in each cluster
#spots_groups$`1`[1:10]
#spots_groups$`2`[1:10]

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

#bsrinf.redP  <- reduceToPathway(bsrinf)  
#bsrinf.redPBP <- reduceToBestPathway(bsrinf.redP) 

bsrsig.redBP <- getLRGeneSignatures(bsrinf.redBP, qval.thres=0.001)
scoresLR <- scoreLRGeneSignatures(bsrdm, bsrsig.redBP, name.by.pathway=FALSE)

simpleHeatmap(scoresLR[1:20,], path=plot_dir, filename=paste0("/",sampleId, "_scoresLR"), column.names=TRUE, height=5, width=9, pointsize=10, hcl.palette="Cividis")

bsrsig.redPBP <- getLRGeneSignatures(bsrinf.redPBP, qval.thres=0.001)
     
scoresPathway <- scoreLRGeneSignatures(bsrdm, bsrsig.redPBP, name.by.pathway=TRUE)
     
#simpleHeatmap(scoresPathway[1:10,], 
#                   path = plot_dir,
#                   filename = paste0("/",sampleId, "_scoresPathway"),
#                   column.names = TRUE, 
#                   width = 9, 
#                   height = 5, 
#                   pointsize = 12,
#                   hcl.palette = "Blue-Red 2"
#                   )
pdf(file.path(plot_dir, paste0(sampleId, "_scoresPathway.pdf")), width=20, height=8)
pheatmap(scoresPathway[1:10,], fontsize_col=7)
dev.off()
