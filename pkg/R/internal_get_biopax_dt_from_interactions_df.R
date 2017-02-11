######################################## internal_get_biopax_dt_from_interactions_df ########################################
internal_get_biopax_dt_from_interactions_df<-
    function(indf_list_element){
        require(stringr)
        
        #component ids, e.g. controller id
        comp_ids<-
            indf_list_element$id
        comp_classes<-
            indf_list_element$class
        #remove id and class columns
        #class is really not necessary, as it is included in the "name" column
        indf_list_element<-
            indf_list_element %>%
            dplyr::select(-id
                          ,-class)
        
        #split up all pertaining elements of the dataframe
        #into names and classes
        #and make biopax-style dataframes from each such element
        indf_list_element_proc<-
            #take part of interactions_df (control, or otherwise)
            1:ncol(indf_list_element) %>%
            lapply(FUN=function(cndex){
                column<-
                    indf_list_element[[cndex]]
                #if all members of the column are NA values
                #return NULL
                if(all(is.na(column)) | all(column %in% "NA")){
                    return(NULL)
                }
                prop_vect<-
                    colnames(indf_list_element)[cndex] %>%
                    #split by underscore separator (if any)
                    strsplit("_") %>%
                    unlist %>%
                    .[-length(.)] 
                
                if(length(prop_vect[!prop_vect %in% internal_banned_props()])>0){
                    #take first letters of the defining names of a column
                    #and add them to the vector of component ids
                    #this would make a unique identifier for a given pair of columns
                    #(name and xref), which belongs to the same group
                    el_ids<-
                        prop_vect[!prop_vect %in% internal_banned_props()] %>% 
                        substring(first = 1
                                  ,last = 1) %>%
                        paste0(collapse="") %>%
                        paste(comp_ids
                              ,.
                              ,sep="_")
                } else {
                    el_ids<-
                        comp_ids
                }
                #ensure char class
                class(column)<-"character"
                #ensure no char NA
                column[column %in% "NA"]<-
                    NA
                #process column, element by element
                #also supply ids from the ids column of the indf_list_element df
                column_proc<-
                    el_ids %>%
                    seq_along %>%
                    lapply(FUN=function(rindex){
                        result<-
                            internal_string_to_df(string=column[rindex]
                                         ,el_id=el_ids[rindex]
                                         #property vect is required to discern 
                                         #xref from entityref-xref combo
                                         ,prop_vect=prop_vect
                            )
                        return(result)
                    }) %>%
                    do.call(rbind.data.frame
                            ,.)
                
                if(nrow(column_proc)>0) {
                    #make auxiliary dataframe
                    #for each parent id and parent class
                    column_proc_aux_df<-
                        comp_ids %>%
                        seq_along %>%
                        lapply(FUN=function(vindex){
                            #filter rows that are pertaining to the same
                            #main parent comp id
                            temp_df<-
                                column_proc %>%
                                filter(grepl(paste0(comp_ids[vindex]
                                                    ,"(?![[:digit:]])")
                                             ,id
                                             ,fixed=FALSE
                                             ,perl=TRUE))
                            if(nrow(temp_df)>0){
                                result<-
                                    internal_make_aux_df_for_dTable(result_dF=temp_df
                                                           ,parent_id=comp_ids[vindex]
                                                           ,parent_class=comp_classes[vindex]
                                                           ,prop_vect=prop_vect
                                    )
                            } else {
                                result<-NULL
                            }
                            return(result)
                        }) %>%
                        do.call(rbind.data.frame
                                ,.)
                    column_proc<-
                        rbind.data.frame(column_proc
                                         ,column_proc_aux_df)
                }
                return(column_proc)
            }) %>%
            do.call(rbind.data.frame
                    ,.) %>% 
            unique 
        
        return(indf_list_element_proc)
        
    }
######################################## internal_get_biopax_dt_from_interactions_df ########################################
