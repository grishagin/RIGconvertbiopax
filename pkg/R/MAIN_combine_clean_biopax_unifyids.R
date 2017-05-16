######################################## MAIN_combine_clean_biopax_unifyids ########################################
MAIN_combine_clean_biopax_unifyids<-
    function(new_biopax_list
             ,pw_df_path="default"
             ,exclude_id_pattern="inxight_pathways"
    ){
        #' @title
        #' MAIN -- Combine and Clean List of BioPAX Objects
        #' @description 
        #' Merges BioPAX objects, unifies component IDs, replaces certain pathway ids with proper \code{inxight_pathway} ones, 
        #' and cleans BioPAX objects from html tags and garbled symbols. 
        #' @details 
        #' Wrapper for \code{unify_biopax_ids} and \code{biopax_from_dt}.
        #' ID unification is conducted by appending a class of a component to its throughout number.
        #' E.g. for a component of class Protein, #11 from the top of the table, ID will be "Protein11".
        #' @param new_biopax_list List of BioPAX objects.
        #' @param pw_df_path Path to the source table with pathway, id, name, and source values.
        #' @param exclude_id_pattern Exclude pathways with such pattern in their IDs from ID unification.
        
        #' @author 
        #' Ivan Grishagin
        
        #establish initial key parameters, if not supplied by user
        if(pw_df_path=="default"){
            pw_df_path<-
                system.file("extdata"
                            ,"pathways_matched_to_sources_current_version.xlsx"
                            ,package="RIGbiopax")
        }
        
        #load pathway spreadsheet
        pw_df<-
            read_excel_astext(path = pw_df_path
                              #,col_types = rep("text",11)
                              ,sheet = 1) %>%
            filter(!is.na(toxdb.Pathway.ID)
                   ,!is.na(biopax.Pathway.ID)) %>%
            dplyr::select(toxdb.Pathway.ID
                          ,toxdb.Pathway.Name
                          ,biopax.Pathway.ID
                          ,toxdb.Pathway.Name
                          ,Source) %>%
            #make up id-source combinations
            #to avoid excluding the same ids that come from different sources
            mutate(bpid_src=
                       paste0(biopax.Pathway.ID
                              ,Source)
                   ,status=
                       "included") %>%
            arrange(Source
                    ,biopax.Pathway.ID
                    ,toxdb.Pathway.ID) %>%
            unique 
        
        
        #find duplicated ids and report them in a file
        pw_df$status[duplicated(pw_df$bpid_src)]<-
            "excluded"
        bp_dupl_ids<-
            pw_df$biopax.Pathway.ID[duplicated(pw_df$bpid_src)]
        
        try(openxlsx::write.xlsx(filter(pw_df
                                        ,biopax.Pathway.ID %in% bp_dupl_ids)
                                 ,paste0(Sys.Date()
                                         ,"_duplicated_ids.xlsx")
                                 ,col.names=TRUE
                                 ,row.names=FALSE))
        write.table(filter(pw_df
                           ,biopax.Pathway.ID %in% bp_dupl_ids)
                    ,paste0(Sys.Date()
                            ,"_duplicated_ids.txt")
                    ,quote = FALSE
                    ,sep = "\t"
                    ,row.names = FALSE
                    ,col.names = TRUE)
        
        #consider only unique ids
        pw_df<-
            pw_df %>%
            filter(status=="included")
        
        combined_biopax<-
            1:length(new_biopax_list) %>% 
            #first, unify ids for each biopax dt
            #and label them with first two letters of biopax source
            lapply(FUN=function(lindex){
                #get ids for a given source
                temp_pw_df<-
                    pw_df %>%
                    filter(Source %in% names(new_biopax_list)[lindex])
                
                missing_ids<-
                    temp_pw_df$biopax.Pathway.ID[!temp_pw_df$biopax.Pathway.ID %in% new_biopax_list[[lindex]]$dt$id]
                
                if(length(missing_ids)>0){
                    print(head(new_biopax_list[[lindex]]$dt))
                    stop("MAIN_combine_clean_biopax_unifyids: in the biopax list element "
                         ,names(new_biopax_list)[lindex]
                         ," not all desired pathways are present.\n"
                         ,"The following ids are not present:\n"
                         ,missing_ids)
                    
                }
                
                #for those ids, replace biopax ids with inxight pathway ids
                temp_unif_df<-
                    new_biopax_list[[lindex]]$dt %>%
                    mutate(id = 
                               mapvalues(id
                                         ,from = temp_pw_df$biopax.Pathway.ID
                                         ,to = temp_pw_df$toxdb.Pathway.ID)
                           ,property_attr_value = 
                               mapvalues(property_attr_value
                                         ,from = temp_pw_df$biopax.Pathway.ID[temp_pw_df$biopax.Pathway.ID %in% property_attr_value]
                                         ,to = temp_pw_df$toxdb.Pathway.ID[temp_pw_df$biopax.Pathway.ID %in% property_attr_value])) 
                temp_unif_df<-
                    temp_unif_df %>%
                    unify_biopax_ids(idtag=substr(tolower(names(new_biopax_list)[lindex])
                                                  ,start = 1
                                                  ,stop = 2)
                                     ,exclude_id_pattern=exclude_id_pattern
                                     ,exclude_class = "Pathway"
                    )
                return(temp_unif_df)
            }) %>% 
            do.call(rbind.data.frame
                    ,.) %>%
            internal_fix_NOTFOUND %>% 
            unique %>% 
            as.data.table %>% 
            #make biopax
            biopax_from_dt %>% 
            
            #remove duplicate biopax components
            remove_duplicate_biopax_components %>% 
            #annotate with gene ids
            add_symbols_entrezids2biopax %>% 
            #unify (re-label) biopax ids
            unify_biopax_ids(exclude_id_pattern=exclude_id_pattern
                             ,exclude_class = "Pathway")
      
        return(combined_biopax)
        
    }
######################################## MAIN_combine_clean_biopax_unifyids ########################################
