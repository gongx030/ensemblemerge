#' Check package compatability and capability for running methods
#'
#' @param package string of package name to check
#' @param version version of package to request install if package is missing or below version number
#' @param language language that package is in
#' @import dplyr
#' @importFrom reticulate import
#'
#' @return returns NULL
#'
.check_package <- function(object){

	browser()

  if(language == "R"){

    out <- tryCatch({
			packageVersion(package) >= version
		},
		error = function(cond){
			if (attributes(cond)$class[1] == 'packageNotFoundError')
				sprintf('%s package %s (>=%s) needs to be installed', language, package, version) %>% message()
			FALSE	
   	}, silent = TRUE) 

	}else if (language == 'Python'){

		is_pip_available <- system('which pip', ignore.stderr = TRUE, ignore.stdout = TRUE) == 0

		if (!is_pip_available){
			message('pip is not available') 
			out <- FALSE	
		}else{
	    out <- tryCatch({
				browser()
				py_module_available(package)
			}, error = function(cond){
				if (grepl('^ModuleNotFoundError', cond$message))
					sprintf('%s package %s (>=%s) needs to be installed', language, package, version) %>% message()
				FALSE
			}, silent = TRUE)
		}	
	}
	 #reticulate::py_run_string(paste("import ", package, "; version = ", package, ".__version__", sep = ""))
#	  import(package)
          #if(!py$version>=version){
          #  message(paste("Later version for package is recommended: ", package))
          #  #message(paste("Please install package version: ", version, " to use this feature", sep = ""))
          #  input <- readline(prompt=paste("Would you like to install package version: ", version, " to use this feature? (Y or N)", sep = ""))
          #  if(input == "y" | input == "Y"){
          #     reticulate::py_install(package, pip = TRUE)
          #  }
          #}
#        },
#        error=function(cond) {
 #           message(paste("Package does not seem to be installed: ", package))
  #          #message(paste("Please install the packages to use to this feature: pip install ", package, "==", version, sep = ""))
   #         input <- readline(prompt=paste("Would you like to install package version: ", version, " to use this feature? (Y or N)", sep = ""))
    #        if(input == "y" | input == "Y"){
     #          reticulate::py_install(sprintf("%s==%s",package, version), pip = TRUE)
			   #system(sprintf("pip install %s==%s", package, version))
#            }
      #  },
#        warning=function(cond) {
 #           message(paste("Package does not seem to be installed: ", package))
  #          message(paste("Please install the packages to use to this feature: pip install ", package, "==", version, sep = ""))
   #         input <- readline(prompt=paste("Would you like to install package version: ", version, " to use this feature? (Y or N)", sep = ""))
#            if(input == "y" | input == "Y"){
#               reticulate::py_install(sprintf("%s==%s",package, version), pip = TRUE)
			   #system(sprintf("pip install %s==%s", package, version))
#            }
#        },
#    )
  return(out)
}
