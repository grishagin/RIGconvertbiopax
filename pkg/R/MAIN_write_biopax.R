######################################## MAIN_write_biopax ########################################
MAIN_write_biopax<-
    function(combined_biopax
             ,what=c("both","onebig","manysmall")
             ,onebigfilename=NULL
             ,pwtoextract_pattern="inxight_pathways"
             ,pwtoextract="default"
             ,bigtarname="default"
             ,smalltarname="default"
    ){
        #' @title
        #' MAIN -- Write BioPAX
        #' @description 
        #' Wrapper for a \code{biopax_from_dt} or \code{splitbiopax_from_dt} functions. 
        #' @details 
        #' Returns either a graph or a dataframe with control-type components.
        #' @param combined_biopax BioPAX object.
        #' @param what How to write: \code{onebig} (all pathways in one file), \code{manysmall} (one pathway per file) or \code{both}.
        #' @param onebigfilename In case one file with all pathways desired, which name should be used for it?
        #' @param pwtoextract_pattern Only pathways with IDs containing the pattern will be extracted.
        #' @param pwtoextract IDs of pathways to extract. Trumps \code{pwtoextract_pattern}.
        #' @param bigtarname Name of the tar.gz archive with one large pathway file.
        #' @param smalltarname Name of the tar.gz archive with individual pathway files.
        
        
        #' @author 
        #' Ivan Grishagin
        
        if(what %in% c("both","onebig")){
            if(onebigfilename=="default"){
                onebigfilename<-
                    paste0(Sys.Date()
                           ,"_full_biopax.owl")
            }
            #single file
            st <- Sys.time()
            new_full_biopax<-
                biopax_from_dt(combined_biopax$dt
                               ,onebigfilename)
            et <- Sys.time()
            
            msg<-
                message("Writing BioPAX object to one large file took "
                        ,round(x = difftime(et,st,units="secs")
                               ,digits = 0)
                        ," seconds.")
            write(msg
                  ,"time_log.txt"
                  ,sep = "\n"
                  ,append=TRUE)
            
            #try to archive the file
            if(bigtarname=="default"){
                bigtarname<-
                    onebigfilename %>% 
                    gsub("\\.owl$"
                         ,""
                         ,.) %>% 
                    paste0("--single_file.tar.gz")
            }
            try(make_targz(input_names = onebigfilename
                           ,output_name = bigtarname))
        }
        
        if(what %in% c("both","manysmall")){
            #individual files
            
            #first, make directory for them and set it as working
            dir.create(path = "./individual_pathways"
                       ,showWarnings = FALSE)
            setwd("./individual_pathways")
            
            #pathways to extract and write
            #if not specified,
            if(pwtoextract=="default"){
                #then search for pathways by pattern
                if(!is.null(pwtoextract_pattern)){
                    pwtoextract<-
                        listPathways(combined_biopax)$id %>%
                        .[grep(pwtoextract_pattern
                               ,.)]
                } else {
                    message("Hey, you didn't specify either pattern or pathway ids to extract!")
                    return()
                }
            }

            st <- Sys.time()
            splitbiopax_list<-
                splitbiopax_from_dt(biopax = combined_biopax
                                    ,write_to_files=TRUE
                                    ,pwtoextract=pwtoextract
                                    ,unifyids=TRUE
                )
            et <- Sys.time()
            msg<-
                message("Writing BioPAX object to multiple files took "
                        ,round(x = difftime(et,st,units="mins")
                               ,digits = 1)
                        ," minutes.")
            write(msg
                  ,"time_log.txt"
                  ,sep = "\n"
                  ,append=TRUE)
            message(msg)
            
            #set working directory back to the original level
            setwd("..")
            
            #try to archive the file
            if(smalltarname=="default"){
                #try to archive
                smalltarname<-
                    paste(Sys.Date()
                          ,"inxight_pathways--individual_files.tar.gz"
                          ,sep="_")
            }
            
            smallfiles<-
                list.files(path = "./individual_pathways"
                           ,pattern="\\.owl$") %>% 
                .[!. %in% onebigfilename]
            try(make_targz(input_names = smallfiles
                           ,output_name = smalltarname))
            
        }
        
    }
######################################## MAIN_write_biopax ########################################
