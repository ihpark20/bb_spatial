
# read tissue positions data (visium)
get_tissue_positions <- function(tissue_positions_list_csv_file) {
    tissue_positions <- read.table(tissue_positions_list_csv_file, sep=",")
    colnames(tissue_positions) <- c("Barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
    rownames(tissue_positions) <- tissue_positions$Barcode    
    tissue_positions <- subset(x=tissue_positions, select = -c(Barcode))
    tissue_positions
}


# https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/spatial
# spot과 spot 간의 distance matrix (unit = micro meter)
# scalefactos 데이터를 이용하여 distance 단위를 micro meter로 변경함
get_spot_distance_matrix <- function(tissue_positions_df, scalefactors_json) {

    spot_distance <- dist(tissue_positions_df[c("pxl_row_in_fullres", "pxl_col_in_fullres")], method="euclidean", diag=TRUE, upper=TRUE)
    spot_distance_matrix <- as.matrix(spot_distance)
    colnames(spot_distance_matrix) <- rownames(tissue_positions_df)
    rownames(spot_distance_matrix) <- rownames(tissue_positions_df)
    scalefactors_data <- jsonlite::fromJSON(scalefactors_json)

    # 65 um (diameter)
    diameter <- 65
    spot_diameter_fullres <- scalefactors_data$spot_diameter_fullres

    # 1 pixel distance ( micron / pixel)
    pixel_distance <-  diameter / spot_diameter_fullres

    # center - center distance (from spot)
    spot_distance_matrix_um <- spot_distance_matrix * pixel_distance
    spot_distance_matrix_um
}




# source spots 에 있는 각 spot에 대하여 target_spots에 가장 가까운 거리
get_min_distance_from_source_spots_to_target_spots <- function(spot_distance_matrix, source_spots, target_spots) {
    sub_matrix <- spot_distance_matrix[source_spots, target_spots]
    min_distance <- apply(sub_matrix[,target_spots], 1, min)
    min_distance
}


# source spots 에 있는 각 spot에 대하여 distance_threshold 이내에 있는 target spots의 개수
get_neighboring_target_spots_count <- function(spot_distance_matrix, source_spots, target_spots, distance_threshold) {
    
    target_spots_number <- apply(spot_distance_matrix[source_spots, target_spots] < distance_threshold, 1, sum)
    target_spots_number
}


