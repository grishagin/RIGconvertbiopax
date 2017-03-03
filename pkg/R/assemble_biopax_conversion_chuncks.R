assemble_biopax_conversion_chuncks<-
    function(bindex
             ,wkdir=NULL
             ,output_dir=NULL){
        #' @title
        #' Assemble Files from BioPAX Conversion Process
        #' @description 
        #' Assembles all files from BioPAX flattening and re-conversion back to BioPAX, 
        #' when it has been carried out a few pathways at a time, resulting in multiple files per source. 
        #' @param bindex Number of the source BioPAX object: \cr
        #' 1-"BioCarta"\cr
        #' 2-"KEGG"\cr
        #' 3-"NCI-Nature"\cr
        #' 4-"NetPath"\cr
        #' 5-"Wiki Pathways"\cr
        #' 6-"Science Signaling"\cr
        #' 7-"Reactome"\cr
        #' 8-"RMC".
        #' @param wkdir Folder with all files.
        #' @param output_dir Where to deposit results.
        

        #' @author 
        #' Ivan Grishagin
        
        if(is.null(wkdir)){
            wkdir<-
                getwd()
        }
        
        if(is.null(output_dir)){
            output_dir<-
                "./merged_results"
        }
        #in case output dir does not exist, create
        dir.create(path = output_dir
                   ,showWarnings = FALSE)
        
        prepareSession(wkdir)
        
        ############# make file lists
        biopax_source_names<-
            c("BioCarta"
              ,"KEGG"
              ,"NCI-Nature"
              ,"NetPath"
              ,"Wiki Pathways"
              ,"Science Signaling"
              ,"Reactome"
              ,"RMC"
            )
        #all files in the folder
        all_files<-
            list.files()
        
        bpname<-
            biopax_source_names[bindex]

        ################# rdata files
        rdata_files<-
            all_files %>% 
            .[grepl(bpname
                    ,.)
              & grepl("\\.RData"
                      ,.)
              & grepl("cleanws"
                      ,.)] 
        sapply(rdata_files
               ,load
               ,envir = environment())
        pieces<-
            grep(ls()
                 ,pattern = paste0(bpname
                                   ,"chunck_biopax")
                 ,value=TRUE)

        #prepare full biopax
        biopax_full<-
            lapply(pieces
                   ,function(biopax_name){
                       biopax_obj<-
                           get(biopax_name)
                       return(biopax_obj$dt)
                   }) %>%
            do.call(rbind.data.frame
                    ,.) %>% 
            unique %>% 
            biopax_from_dt
        #and assign it to properly named variable
        #for future convenience
        biopax_name<-
            paste0(gsub(" |-"
                        ,"_"
                        ,tolower(bpname))
                   ,"_biopax")
        assign(biopax_name
               ,biopax_full)
        
        rdata_outfile<-
            paste(Sys.Date()
                  ,bpname
                  ,"biopax_object_workspace.RData"
                  ,sep="_")
        
        ################# comparison files
        comparison_df<-
            all_files %>% 
            .[grepl(bpname
                    ,.)
              & grepl("comparison"
                      ,.)
              & grepl("\\.txt"
                      ,.)] %>% 
            lapply(read.delim
                   ,quote=""
                   ,check.names = FALSE) %>% 
            do.call(rbind.data.frame
                    ,.) %>% 
            mutate(Source=bpname)
        
        comparison_outfile<-
            paste(Sys.Date()
                  ,bpname
                  ,"comparisons_new2orig.txt"
                  ,sep="_")
        
        ################# interactions
        interaction_df<-
            all_files %>% 
            .[grepl(bpname
                    ,.)
              & grepl("interaction"
                      ,.)
              & grepl("\\.txt"
                      ,.)] %>% 
            lapply(read.delim
                   ,quote="") %>% 
            do.call(rbind.data.frame
                    ,.)
        
        interaction_outfile<-
            paste(Sys.Date()
                  ,bpname
                  ,"interactions.txt"
                  ,sep="_")
        
        
        ################# write to files
        write.table(comparison_df
                    ,file.path(output_dir
                               ,comparison_outfile)
                    ,quote = FALSE
                    ,sep = "\t"
                    ,row.names = FALSE
                    ,col.names = TRUE)
        
        write.table(interaction_df
                    ,file.path(output_dir
                               ,interaction_outfile)
                    ,quote = FALSE
                    ,sep = "\t"
                    ,row.names = FALSE
                    ,col.names = TRUE)
        #gather all appropriate variables
        all_vars<-
            ls(envir = environment())
        #and remove them
        rm(list=all_vars[!all_vars %in% c(biopax_name
                                          ,"rdata_outfile")]
           ,envir = environment())
        
        #save image
        save(list = ls(all.names = TRUE)
             ,file = 
                 file.path(output_dir
                           ,rdata_outfile)
             ,envir = environment())
        return()
    }
