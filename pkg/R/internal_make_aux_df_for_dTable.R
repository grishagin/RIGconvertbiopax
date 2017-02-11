######################################## internal_make_aux_df_for_dTable ########################################
internal_make_aux_df_for_dTable<-
    function(result_dF
             ,parent_id
             ,parent_class
             ,prop_vect){
        #function basically makes auxiliary table 
        #that establishes all proper references 
        #between the very top level component id 
        #and all ids for each sub-component 
        #by establishing connection between properties
        #and the numbers attached to each id 
        prop_vect<-
            prop_vect[!prop_vect %in% internal_banned_props()]
        
        #all ids that have not been referenced
        nonrefids<-
            result_dF$id[!result_dF$id %in% result_dF$property_attr_value &
                             !result_dF$id %in% parent_id] %>%
            unique
        
        if(length(nonrefids)<1){
            return(NULL)
        }
        #make a df from these ids by splitting them into "main" id and markers
        nonrefids_df<-
            idvect_to_df_by_index(nonrefids)
        
        if(is.null(nonrefids_df)){
            ids_list<-
                list(parent_id)
        } else {
            #add properties based on level
            ids_list<-
                prop_vect[-length(prop_vect)] %>%
                seq_along %>%
                lapply(FUN=function(pindex){
                    #paste first letters of properties
                    alpha_vect<-
                        prop_vect[1:pindex] %>%
                        substring(first = 1
                                  ,last = 1) %>%
                        paste(collapse="")
                    #paste level numbers
                    num_vect<-
                        nonrefids_df %>%
                        dplyr::select(1:pindex) %>%
                        apply(MARGIN = 1
                              ,paste
                              ,collapse="-")
                    #paste all those together, add parent id, and return
                    result<-
                        paste(alpha_vect
                              ,num_vect
                              ,sep="-") %>%
                        paste(parent_id
                              ,.
                              ,sep="_")
                    return(result)
                }) %>%
                #join made-up ids with parent id
                #parent id goes first
                c(parent_id
                  ,.)
        }
        
        #if there are any made-up ids, cook classes and pavs differently
        #otherwise the result is undesirable 
        if(length(ids_list)>1){
            #list of classes
            # classes_list<-
            #     list(parent_class
            #          ,result_dF$class[match(nonrefids
            #                                 ,result_dF$id)])
            if(prop_vect[length(prop_vect)]=="component"){
                final_class<-
                    "Complex"
            } else {
                final_class<-
                    result_dF$class[match(nonrefids
                                          ,result_dF$id)]
            }
            classes_list<-
                list(parent_class
                     ,final_class)
            #list of pavs
            pav_list<-
                #made-up ids plus non-referenced ids
                ids_list[-1] %>%
                c(list(nonrefids))
        } else {
            classes_list<-
                parent_class
            pav_list<-
                #just non-referenced ids
                list(nonrefids)
        }
        
        #now go through the lists/vectors of elements
        #and make biopax-style data tables from them
        #by joining them accordingly
        dTable<-
            prop_vect %>%
            seq_along %>%
            lapply(FUN=function(pindex){
                #if there are more than two elements in the 
                #prop_vector, i.e. it's longer than right->component, etc.
                #class could not be determined (at least at this point)
                #but this function is prepared to deal with these cases -- just in case
                #by using a "NOT FOUND" stub
                if(pindex %in% c(1,length(prop_vect))){
                    temp_class<-
                        classes_list[[pindex]]
                } else {
                    temp_class<-
                        "CLASS_NOT_FOUND"
                }
                #make the appropriate data table
                dT<-
                    data.frame(class=temp_class
                               ,id=ids_list[[pindex]]
                               ,property=prop_vect[pindex]
                               ,property_attr="rdf:resource"
                               ,property_attr_value=pav_list[[pindex]]
                               ,property_value=""
                    ) 
                return(dT)
            }) %>%
            do.call(rbind.data.frame
                    ,.) %>%
            unique
        
        return(dTable)
        
    }
######################################## internal_make_aux_df_for_dTable ########################################
