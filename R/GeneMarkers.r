setClass(
	'PanglaoDBGeneMarkers',
	representation(
		celltype = 'data.frame',
		genome = 'character',
		level = 'character'
	),
	contains = 'BaseGeneMarkers',
	prototype(
		name = 'PanglaoDBGeneMarkers',
		check_dependencies = TRUE,
		genome = 'hg19',
		level = 'Level4',
		dependences = list(
			new('RPackage', package_name = 'scMRMA', package_version = '1.0')
		)
	),
	validity = function(object){
		msg <- NULL
		if (!object@genome %in% c('hg19', 'mm10'))
			msg <- sprintf('unknown genome: %s', object@genome)
		return(msg)
	}
)

#' @importFrom methods callNextMethod
#'
setMethod('initialize', 'PanglaoDBGeneMarkers', function(.Object, check_dependencies = TRUE, ...){
	.Object <- callNextMethod(.Object, check_dependencies = check_dependencies, ...)
	if (.Object@genome == 'hg19'){
		.Object@celltype <- local({get(load(system.file("data", "Human_PanglaoDB.Rdata", package = "scMRMA")))})
	}else if (.Object@genome == 'mm10'){
		.Object@celltype <- local({get(load(system.file("data", "Mouse_PanglaoDB.Rdata", package = "scMRMA")))})
	}
	.Object
})
