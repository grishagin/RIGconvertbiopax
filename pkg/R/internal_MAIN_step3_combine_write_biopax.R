internal_MAIN_step3_combine_write_biopax<-
    function(file_dir=NULL
             ,output_dir=NULL){
        
        if(is.null(file_dir)){
            file_dir<-
                getwd()
        }
        
        prepareSession(file_dir)
        
        if(is.null(output_dir)){
            output_dir<-
                "./combined_biopax"
            dir.create(path = output_dir
                       ,showWarnings = FALSE)
        }

        setwd(output_dir)
        
        all_files<-
            list.files(path = file_dir
                       ,full.names = TRUE)
        
        ############# src names
        #PAY ATTENTION TO ORDER!
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
        biopax_obj_names<-
            c("biocarta_biopax"
              ,"kegg_biopax"
              ,"nci_nature_biopax"
              ,"netpath_biopax"
              ,"wiki_pathways_biopax" 
              ,"science_signaling_biopax"
              ,"reactome_biopax"
              ,"rmc_biopax"
            )
        
        
        ################# rdata files
        #present sources
        present_biopax_logi<-
            biopax_obj_names %in% ls(envir = environment())
        #redefine source names 
        biopax_source_names<-
            biopax_source_names[present_biopax_logi]
        biopax_obj_names<-
            biopax_obj_names[present_biopax_logi]
            
        rdata_files<-
            all_files %>% 
            .[grepl("biopax_object_workspace"
                    ,.)] 
        invisible(sapply(rdata_files
                         ,load
                         ,envir = environment()))
        
        #prepare biopax objects' list
        biopax_obj_list<-
            biopax_obj_names %>% 
            lapply(get
                   ,envir=environment())
        names(biopax_obj_list)<-
            biopax_source_names
        ################# comparison files
        comparison_df<-
            all_files %>% 
            .[grepl("comparison"
                    ,.)] %>% 
            lapply(read.delim
                   ,quote=""
                   ,check.names = FALSE) %>% 
            do.call(rbind.data.frame
                    ,.) 
        
        comparison_outfile<-
            paste(Sys.Date()
                  ,"ALL_comparisons_new2orig.txt"
                  ,sep="_")
        
        ################# write to files
        write.table(comparison_df
                    ,file.path(comparison_outfile)
                    ,quote = FALSE
                    ,sep = "\t"
                    ,row.names = FALSE
                    ,col.names = TRUE)
        
        
        
        combined_biopax<-
            MAIN_combine_clean_biopax_unifyids(new_biopax_list=biopax_obj_list)
        
        rm(list=ls()[!ls() %in% "combined_biopax"]
           ,envir = environment())
        
        save.image(paste0(Sys.Date()
                          ,"_combined_biopax.RData"))
        
        #write biopax to file
        onebigfilename<-
            paste0(Sys.Date()
                   ,"_all_inxight_pathways.owl")
        MAIN_write_biopax(combined_biopax=combined_biopax
                          ,what="both"
                          ,pwtoextract_pattern="inxight_pathways"
                          ,onebigfilename=onebigfilename)
    }

