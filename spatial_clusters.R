library(config)
library(BulkSignalR)
library(igraph)
library(Seurat)
library(ggpubr)

# Load config data
config <- config::get(file="config.yaml")
data_dir <- config$data
proxy_dir <- config$proxy
output_dir <- config$output
plot_dir <- file.path(output_dir, "immuneSpatial", "figures")

source("utils.R")

sampleId <- "YS15-33091"
seurat_processed <- file.path(proxy_dir, paste0(sampleId, "_seurat_processed.RDS"))
seurat_obj <- readRDS(seurat_processed)

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

kk1 <- FindMarkers(seurat_obj, ident.1 = tumor_near_spots, ident.2 = tumor_dist_spots, only.pos=TRUE)
kk2 <- FindMarkers(seurat_obj, ident.1 = tumor_periphery_spots, ident.2 = tumor_near_spots, only.pos=TRUE)
kk3 <- FindMarkers(seurat_obj, ident.1 = tumor_core_spots, ident.2 = tumor_periphery_spots, only.pos=TRUE)

kk1_only <- setdiff(rownames(kk1), rownames(kk2))
kk1_only <- setdiff(kk1_only, rownames(kk3))

kk1[kk1_only,]

tx <- GetAssay(seurat_obj)[kk1_only,]
tx <- data.frame(dx)

t_cell_correlation <- apply(tx, 1, function(x) cor(as.numeric(x), seurat_obj$T.cells_ratio))
b_cell_correlation <- apply(tx, 1, function(x) cor(as.numeric(x), seurat_obj$B.cells_ratio))

#distances_from_tumor <- get_min_distance_from_source_spots_to_target_spots(spot_dist_matrix, spots_names_selected, tumor_spots) 

#t_cell_ratio_thr <- 0.3
#spots_selected <- seurat_obj$T.cells_ratio > t_cell_ratio_thr
#spots_names_selected <- names(spots_selected[spots_selected==T])

#table(seurat_obj$patho[spots_names_selected])
#SpatialDimPlot(seurat_obj, cells.highlight=spots_names_selected, label=F) +  theme(legend.position = "none")
#dev.off()



near_spots <- names(distances_from_tumor[distances_from_tumor < 150])
far_spots <- names(distances_from_tumor[distances_from_tumor >= 150])

kk <- FindMarkers(seurat_obj, ident.1 = near_spots, ident.2 = far_spots)
kk <- FindMarkers(seurat_obj, ident.1 = near_spots, ident.2 = far_spots, only.pos=TRUE)


# defining tumor core and tumor boundary
# tumor core: distance from non-tumor spots > 150 (?)
# tumor boundary: distance from non-tumor spots < 150 (?)
tumor_core <- names(distances_from_tumor[distances_from_tumor < 50])





# distanace matrix for the selected spots
spot_dist_selected <- spot_dist_matrix[spots_names_selected, spots_names_selected]


# sptailly clustered spots (looking for neighbors)
# make igraph

distance_thr <- 160
g <- graph.adjacency(spot_dist_selected < distance_thr, mode="undirected", weighted=T, diag=F)
summary(g)

results_components <- components(g)
unique_components <- unique(results_components$membership)

membership <- results_components$membership
patho_anno <- seurat_obj$patho[spots_names_selected]

min_spots_in_component <- 5


# define your interest gruop
tumor_annotations <- c("DCIS", "invasive carcinoma", "tumor", "cancer")
tumor_spots <- names(seurat_obj$patho[seurat_obj$patho %in% tumor_annotations])

for (component_number in unique_components) {
    component_spots <- names(membership[membership == component_number])
    if (length(component_spots) < min_spots_in_component) {
        next
    }
    message(component_number, ",", length(component_spots))
    #print(seurat_obj$T.cells[component_spots])

    component_dist_to_tumor_spots <- get_min_distance_from_source_spots_to_target_spots(spot_dist_matrix, component_spots, tumor_spots) 
    print(mean(component_dist_to_tumor_spots))
    print("--")
}

t_spots <- seurat_obj[,spots_names_selected]

t_spots$membership <- membership

selected_components <- results_components$csize >= 5
selected_components <- which(selected_components)

t_spots_connected <- t_spots[,t_spots$membership %in% selected_components]

Idents(t_spots_connected) <- "membership"
kk <- FindAllMarkers(t_spots_connected)


#vis_data_selected <- vis_data[,tissue_positions_df_in_tissue_tme$Barcode]
#Idents(vis_data_selected) <- tissue_positions_df_in_tissue_tme$DIST_ZONE

#deg_zone <- FindAllMarkers(vis_data_selected, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)

kk %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10

pdf(file.path(plot_dir, paste0(sampleId, "_tcell_spots_heatmap.pdf")), width=10, height=15)
DoHeatmap(t_spots_connected, features = top10$gene, label=TRUE, size=2, angle=30) 
dev.off()




#t_cells_df <- data.frame(T.cells=seurat_obj$T.cells, T.cells_ratio=seurat_obj$T.cells_ratio)
#pdf(file.path(plot_dir, paste0(sampleId, "_T.cells_ratio.pdf")))
#ggdensity(t_cells_df, x="T.cells_ratio")
#dev.off()
