internal_MAIN_step4_plot_pathways<-
    function(file_dir=NULL
             ,output_dir=NULL
             ,pw_id_pattern="inxight_pathways"){
        
        #pipeline step 4
        #produce SVG files
        
        #prepare directories if not supplied
        if(is.null(file_dir)){
            file_dir<-
                getwd()
        }
        
        if(is.null(output_dir)){
            output_dir<-
                file.path(file_dir
                          ,"combined_biopax")
            dir.create(path = output_dir
                       ,showWarnings = FALSE)
        }
        
        prepareSession(output_dir)
        
        ################# load the rdata file
        rdata_file<-
            list.files(path = file_dir
                       ,full.names = TRUE
                       ,pattern = "RData$")
      
        invisible(load(rdata_file
                       ,envir = environment()))
        
        ################# form a list pathways
        if(!"combined_biopax" %in% ls(envir = environment())){
            stop("internal_MAIN_step4_all_biopax_assemble:
                 can't seem to find a 'combined_biopax' object! Aborting...")
        }
        
        pw_ids<-
            listPathways(combined_biopax
                         ,) %>% 
            as.data.table %>% 
            .[grepl(pw_id_pattern
                    ,id)]$id
        
        ################# plot those pathways
        RIGplotbiopax:::plot_pathways(biopax = combined_biopax
                                      ,pw_ids = pw_ids)
        
        return(NULL)
    }


