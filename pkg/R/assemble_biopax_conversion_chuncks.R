assemble_biopax_conversion_chuncks<-
    function(folder
             ,bindex){
        #' @title
        #' Assemble Files from BioPAX Conversion Process
        #' @description 
        #' Assembles all files from BioPAX flattening and re-conversion back to BioPAX, 
        #' when it has been carried out a few pathways at a time, resulting in multiple files per source. 
        #' @param folder Folder with all files.
        #' @param orig_biopax Number of the source BioPAX object: \cr
        #' 1-"BioCarta"\cr
        #' 2-"KEGG"\cr
        #' 3-"NCI-Nature"\cr
        #' 4-"NetPath"\cr
        #' 5-"Wiki Pathways"\cr
        #' 6-"Science Signaling"\cr
        #' 7-"Reactome"\cr
        #' 8-"RMC".

        #' @author 
        #' Ivan Grishagin
        
        dir.create(path = "./merged_results"
                   ,showWarnings = FALSE)
        
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
               ,envir = .GlobalEnv)
        pieces<-
            grep(ls()
                 ,pattern = paste0(bpname
                                   ,"chunck_biopax")
                 ,value=TRUE)
        biopax_full<-
            lapply(pieces
                   ,function(biopax_name){
                       biopax_obj<-
                           get(biop_name)
                       return(biopax_obj$dt)
                   }) %>%
            do.call(rbind.data.frame
                    ,.) %>% 
            biopax_from_dt
        
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
                      ,.)] %>% 
            lapply(read.delim
                   ,quote="") %>% 
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
                    ,file.path("./merged_results"
                               ,comparison_outfile)
                    ,quote = FALSE
                    ,sep = "\t"
                    ,row.names = FALSE
                    ,col.names = TRUE)
        
        write.table(interaction_df
                    ,file.path("./merged_results"
                               ,interaction_outfile)
                    ,quote = FALSE
                    ,sep = "\t"
                    ,row.names = FALSE
                    ,col.names = TRUE)
        
        rm(list=ls()[!ls() %in% c("biopax_full"
                                  ,rdata_outfile)])
        
        #save image
        save.image(file.path("./merged_results"
                             ,rdata_outfile))
           
    }
