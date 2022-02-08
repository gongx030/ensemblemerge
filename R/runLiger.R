#' Run Liger merging
#'
#' @param params a LigerParams object
#' @param data a Seurat object
#'
#' @return returns a Seurat object of the integrated data
#'
run_Liger <- function(params, data){

  data <- data %>% 
		SeuratWrappers::RunOptimizeALS(
			k = params@npcs, 
			lambda = params@lambda, 
			split.by = params@batch, 
			reduction = params@name, 
			reduction.name = params@name, 
			reduction.key = params@name
		) %>%
 		 	SeuratWrappers::RunQuantileNorm(
				split.by = params@batch, 
				reduction = params@name, 
				reduction.name = params@name, 
				reduction.key = params@name
			)
	data
}
