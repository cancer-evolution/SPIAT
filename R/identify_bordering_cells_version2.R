#' identify_bordering_cells_version2
#'
#' @description Identifies the cells bordering a group of cells of a particular phenotype
#'
#' @param sce_object SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param reference_marker Cells positive for this marker will be used as reference
#' @param radius Number specifying the radius for computing clusters. 
#' @param member Number specifying the smallest number of members to be kept in a cluster. Larger number, fewer clusters left.
#' @param alpha Number specifying the ahull parameter. Larger number, more points included in the ahull.
#' @import SingleCellExperiment
#' @export

#sce_object <- sce_ovarian_panimmune1
#reference_marker <- "WT1"
#member <- 100
#alpha <- 50 #larger the more lenient it is
#radius <- 60 #larger the more cells are considered

# colData() is in package 'SummarizedExperiment' but imported by SingleCellExperiment

identify_bordering_cells_version2 <- function(sce_object, reference_marker, radius, alpha, member) {
  formatted_data <- data.frame(colData(sce_object))
  reference_cells <- formatted_data[formatted_data$Phenotype == reference_marker, ]
  clusters <- identify_cell_clusters(sce_object, reference_marker, radius)
  clusters <- clusters[,c("Cell.X.Position", "Cell.Y.Position","Cluster")]
  
  #CHECK
  if (nrow(reference_cells) == 0){
    stop("There are no reference_cells in the dataset for the specified marker")
  }

  table_cluster <- table(clusters$Cluster)
  
  border_cells = c()
  for (i in names(table_cluster[table_cluster > member])){
    if (i != "Cluster_NA"){
      cluster = clusters[clusters$Cluster == i,]
      cluster = unique(cluster)
      ahull_cluster = ahull(cluster$Cell.X.Position, 
                            cluster$Cell.Y.Position, alpha = alpha)
      borders_cluster = cbind(data.frame(ahull_cluster$ashape.obj$edges)$x1,data.frame(ahull_cluster$ashape.obj$edges)$y1)
      border_cells = cbind(border_cells, t(borders_cluster))
    }
  }
  border_cells = data.frame(t(border_cells))
  #CHECK
  if (nrow(border_cells) == 0){
    stop("There are no border cells found")
  }
  plot(border_cells, cex = 0.3)

}
