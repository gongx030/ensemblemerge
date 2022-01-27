run_Seurat <- function(params, data){

	cn <- colnames(data)

  batch_list <- SplitObject(data, split.by = params@batch)

  for(i in 1:length(batch_list)) {
    if(params@norm_data){
      batch_list[[i]] <- NormalizeData(
				batch_list[[i]],
				normalization.method = params@norm_method, 
				scale.factor = params@scale_factor
			)
    }
		batch_list[[i]] <- FindVariableFeatures(
			batch_list[[i]], 
			selection.method = params@selection.method, 
			nfeatures = params@numHVG
		)
  }

  cell_anchors <- FindIntegrationAnchors(object.list = batch_list, dims = 1:params@npcs)
  data <- IntegrateData(anchorset = cell_anchors, dims = 1:params@npcs, k.weight = params@k.weight)

  if(params@regressUMI && params@scaling) {
    data <- ScaleData(object = data , vars.to.regress = params@vars.to.regress)
  } else if(params@scaling) { # default option
    data <- ScaleData(object = data)
  }

  data <- RunPCA(data, npcs = params@npcs, reduction.name = params@name)

	data

}

