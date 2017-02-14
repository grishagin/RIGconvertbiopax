update_RIGconvertbiopax <-
    function () 
    {
        #' @title
        #' Update RIGconvertbiopax from GitHub
        #' @description 
        #' Updates the package. It may be necessary to restart the R session after the update.
        
        #' @author 
        #' Ivan Grishagin
        
        
        #unloadNamespace(ns = "RIGconvertbiopax")
        devtools::unload(pkg = inst("RIGconvertbiopax"))
        devtools::install_github("grishagin/RIGconvertbiopax", subdir = "pkg")
    }