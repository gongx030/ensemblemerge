#' Check package compatability and capability for running methods
#'
#'
#' @param package string of package name to check
#' @param version version of package to request install if package is missing or below version number
#' @param language language that package is in
#'
#' @return returns NULL
#' @export
checkPackage <- function(package, version, language = "R") {
  if(language == "R"){
    out <- tryCatch(
        {
          if(!packageVersion(package)>=version){
            message(paste("Later version for package is recommended: ", package))
            message(paste("Please install package version: ", version, " to use this feature", sep = ""))
          }
        },
        error=function(cond) {
            message(paste("Package does not seem to be installed: ", package))
            message(paste("Please install the packages to use to this feature: install.packages('", package, "')", sep = ""))
        },
        warning=function(cond) {
            message(paste("Package does not seem to be installed: ", package))
            message(paste("Please install the packages to use to this feature: install.packages('", package, "')", sep = ""))
        },
        finally={
        }
    )
  }
  else if(language == "Python"){
    out <- tryCatch(
        {
	  reticulate::py_run_string(paste("import ", package, "; version = ", package, ".__version__", sep = ""))
          if(!py$version>=version){
            message(paste("Later version for package is recommended: ", package))
            message(paste("Please install package version: ", version, " to use this feature", sep = ""))
          }
        },
        error=function(cond) {
            message(paste("Package does not seem to be installed: ", package))
            message(paste("Please install the packages to use to this feature: pip install ", package, "==", version, sep = ""))
        },
        warning=function(cond) {
            message(paste("Package does not seem to be installed: ", package))
            message(paste("Please install the packages to use to this feature: pip install ", package, "==", version, sep = ""))
        },
        finally={
        }
    )
  }
  else{
    message("Language options are either 'R' or 'Python")
  }
  return(out)
}