internal_MAIN_step1_all_biopax_split_interaction_tables<-
    function(groupnum = groupnum 
             ,ngroups = ngroups
             ,bindex = bindex
             ,work_dir = work_dir){
        
        ########### work dir
        prepareSession(work_dir)
        
        ########### workspace with biopax
        biopax_ws<-
            system.file("extdata"
                        ,"biopax_objects_current_version.RData"
                        ,package="RIGbiopax")
        load(biopax_ws)
        pwid_to_convert = "all"
        
        ########### pw df
        
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
        biopax_obj<-
            list(biocarta_biopax
                 ,kegg_biopax
                 ,nci_biopax
                 ,netpath_biopax
                 ,wiki_biopax
                 ,scisig_biopax
                 ,reactome_biopax
                 ,rmc_biopax
            )
        
        ########################## convert biopax to flat file
        temp_inter_df<-
            MAIN_flatten_biopax(pwid_to_convert = pwid_to_convert
                                ,work_dir=work_dir
                                ,pw_df_path="default"
                                ,biopax=biopax_obj[[bindex]]
                                ,source_name=biopax_source_names[bindex]
                                ,groupnum=groupnum
                                ,ngroups=ngroups)
        
        #assign to a variable that has source name and group number in it
        df_var_name<-
            paste0(biopax_source_names[bindex]
                   ,"chunck_df"
                   ,groupnum)
        assign(df_var_name
               ,temp_inter_df
               ,envir = environment())
        #save image
        save(list = ls(all.names = TRUE)
             ,file = 
                 paste0("df_"
                        ,biopax_source_names[bindex]
                        ,groupnum
                        ,".RData")
             ,envir = environment())

        temp_new_biopax<-
            temp_inter_df %>% 
            interactions_df_to_biopax_dt %>%
            biopax_from_dt(filename=NULL)
        
        #assign to a variable that has source name and group number in it
        bp_var_name<-
            paste0(biopax_source_names[bindex]
                   ,"chunck_biopax"
                   ,groupnum)
        assign(bp_var_name
               ,temp_new_biopax
               ,envir = environment())
        
        save(list = ls(all.names = TRUE)
             ,file = 
                 paste0("biopax_df_"
                        ,biopax_source_names[bindex]
                        ,groupnum
                        ,".RData")
             ,envir = environment())
        
        try(MAIN_biopax_comparison(pwid_to_compare=pwid_to_convert
                                   ,pwid_to_plot="none"
                                   ,new_biopax=temp_new_biopax
                                   ,orig_biopax=biopax_obj[[bindex]]
                                   ,verbose = FALSE
                                   ,pw_df_path="default"
                                   ,source_name=biopax_source_names[bindex]
                                   ,groupnum=groupnum
                                   ,ngroups=ngroups))
        
        clean_save_name<-
            paste0("cleanws_biopax_df_"
                   ,biopax_source_names[bindex]
                   ,groupnum
                   ,".RData")
        #do one final, clean save w/o extra variables
        rm(list=ls()[!ls() %in% c(df_var_name
                                  ,bp_var_name
                                  ,"clean_save_name")]
           ,envir = environment())

        save(list = ls(all.names = TRUE)
             ,file = 
                 clean_save_name
             ,envir = environment())
    }
      