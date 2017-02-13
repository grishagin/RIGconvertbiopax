######################################## MAIN_write_biopax ########################################
MAIN_write_biopax<-
  function(combined_biopax
           ,what=c("both","onebig","manysmall")
           ,onebigfilename=NULL
           ,pwtoextract_pattern="inxight_pathways"
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
 
    #' @author 
    #' Ivan Grishagin

    if(what %in% c("both","onebig")){
      if(is.null(onebigfilename)){
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
        paste0("The whole thing took "
               ,round(x = difftime(et 
                                   ,st
                                   ,units = "mins")
                      ,digits = 0)
               ," minutes.")
      write(msg
            ,"time_log.txt"
            ,sep = "\n"
            ,append=TRUE)
      message(msg)
    }
    
    if(what %in% c("both","manysmall")){
      #individual files
      
      #pathways to extract and write
      pwtoextract<-
        listPathways(combined_biopax)$id %>%
        .[grep(pwtoextract_pattern
               ,.)]
      
      st <- Sys.time()
      splitbiopax_list<-
        splitbiopax_from_dt(biopax = combined_biopax
                            ,write_to_files=TRUE
                            ,pwtoextract=pwtoextract
                            ,unifyids=TRUE
        )
      et <- Sys.time()
      msg<-
        paste0("The whole thing took "
               ,round(x = difftime(et 
                                   ,st
                                   ,units = "mins")
                      ,digits = 0)
               ," minutes.")
      write(msg
            ,"time_log.txt"
            ,sep = "\n"
            ,append=TRUE)
      message(msg)
    }
  }
######################################## MAIN_write_biopax ########################################
