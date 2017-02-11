######################################## interactions_df_to_biopax_dt ########################################
interactions_df_to_biopax_dt<-
    function(interactions_df){
        st<-Sys.time()
        message("Converting a dataframe with interactions into a biopax data table...")
        
        #prepare auxiliary dataframe with pw names, ids, etc.
        aux_df<-
            interactions_df %>%
            internal_make_main_aux_df
        
        #prepare interactions_df for processing and extraction of components
        interactions_df_list<-
            list() 
        #controllers...
        interactions_df_list$controller<-
            interactions_df %>%
            filter(!is.na(controller_id)) %>%
            .[,grep("controller_"
                    ,colnames(.))] 
        
        #... controlled...
        interactions_df_list$controlled<-
            interactions_df %>% 
            filter(!is.na(controlled_id)) %>%
            .[,grep("controlled_"
                    ,colnames(.))] 

        #... and other non-control components 
        #ATTENTION: non_ctrl_component column must be logical
        if (class(interactions_df$non_ctrl_component)!="logical"){
            interactions_df$non_ctrl_component<-
                as.logical(interactions_df$non_ctrl_component)
        }
        interactions_df_list$non_ctrl_component<-
            interactions_df %>%
            filter(non_ctrl_component) %>%
            .[,grep("non_ctrl_component_|^component_"
                    ,colnames(.))] 
        
        #strip column name identifiers to make them uniform regardless of exact case
        for(index in 1:length(interactions_df_list)){
            colnames(interactions_df_list[[index]])<-
                colnames(interactions_df_list[[index]]) %>%
                gsub("controller_|controlled_|non_ctrl_component_|component_"
                     ,""
                     ,.) %>%
                gsub("complex"
                     ,"component_"
                     ,.)
        }
        
        #get dt of all components
        comp_df<-
            lapply(interactions_df_list
                   ,internal_get_biopax_dt_from_interactions_df) %>%
            do.call(rbind.data.frame
                    ,.) 
        
        #assemble the final dt
        dTable<-
            rbind.data.frame(aux_df
                             ,comp_df)
        
        #remove duplicates
        #remove rows with NA values
        #replace placeholder classes
        dTable<-
            dTable %>%
            unique %>%
            .[complete.cases(.),] %>% 
            #replace PLACEHOLDERs 
            #with corresponding classes and properties
            internal_fix_placeholders
       

        et<-Sys.time()
        
        msg<-
            paste0("Conversion of a flat file with interactions took "
                   ,round(difftime(et
                                   ,st
                                   ,units = "mins")
                          ,1)
                   ," minutes.")
        message(msg)
        write(msg
              ,"log.txt"
              ,sep = "\n"
              ,append=TRUE)
        
        return(dTable)
    }
######################################## interactions_df_to_biopax_dt ########################################
