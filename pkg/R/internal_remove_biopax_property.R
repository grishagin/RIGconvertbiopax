internal_remove_biopax_property<-
    function(biopax
             ,property_to_remove){
        
        #remove the rest of the memberPhysicalEntity entities
        #entirely by connecting them to their targets directly
        #1. get instances with property
        #2. of those, take only refids
        biopax$dt$property_attr_value<-
            biopax$dt$property_attr_value %>% 
            striphash
            
        prop_df<-
            biopax$dt %>% 
            filter(property==property_to_remove
                   ,property_attr_value %in% id) %>% 
            dplyr::select(id
                          ,refid=property_attr_value)

        
        #3. prepare a list of ids to remove later
        #keep all those ids that are referenced by the referenced ids
        #i.e. if memberPhysicalEntity references a protein -- 
        #all those ids, which that protein references, need to be kept
        ids_to_keep<-
            getReferencedIDs(biopax=biopax
                             ,id=prop_df$refid
                             ,recursive=TRUE)  %>% 
            c(prop_df$refid)
        ids_to_remove<-
            getReferencedIDs(biopax=biopax
                             ,id=prop_df$id
                             ,recursive=TRUE) %>% 
            c(prop_df$id) %>% 
            .[!. %in% ids_to_keep]
        
        #4. merge refids id-wise and then 
        #for each instance referring to the desired property instance
        #replace refid with merged ids
        prop_df<-
            prop_df %>% 
            #but first also add original ids to be able to keep the name properties
            rbind.data.frame(data.frame(id=unique(.$id)
                                        ,refid=unique(.$id))) %>% 
            merge_cols_shorten_df(colKey="id")
        
        biopax$dt$property_attr_value<-
            mapvalues(biopax$dt$property_attr_value
                      ,from = prop_df$id
                      ,to = prop_df$refid) 
        
        #5. remove all instances of undesired property
        #(except those with non-zero value, i.e. name, etc.)
        #and ids they are referring to
        #and expand
        biopax$dt<-
            biopax$dt %>% 
            filter(!(id %in% ids_to_remove
                     & property_value == "")) %>% 
            split_cols_lengthen_df(colsToSplit = "property_attr_value") %>% 
            as.data.table
        
        if("memberPhysicalEntity" %in% biopax$dt$property){
            biopax<-
                biopax %>% 
                internal_remove_biopax_property(property_to_remove="memberPhysicalEntity")
        }
        
        return(biopax)
    }