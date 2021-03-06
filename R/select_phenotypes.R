#' select_phenotypes
#'
#' @description Select phenotypes to keep or exclude in the analysis
#' @param sce_object SingleCellExperiment object generated by format_image_to_sce
#' @param keep TRUE if vector of phenotypes are the cells that are going to be kept,
#' FALSE if they are to be removed
#' @param phentoypes Vector of phenotypes of keep or exclude
#' @export

select_phenotypes <- function(sce_object, keep=TRUE, phenotypes = NULL){
  if(keep){
    sce_object2 <- sce_object[,sce_object$Phenotype %in% phenotypes]
  } else{
    sce_object2 <- sce_object[,!(sce_object$Phenotype %in% phenotypes)]
  }
  return(sce_object2)
}
