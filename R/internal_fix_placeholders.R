######################################## internal_fix_placeholders ########################################
internal_fix_placeholders<-
    function(dF
             ,notfound="NOTFOUND"){
        require(qdap)
        #called by internal_get_biopax_dt_from_interactions_df
        #dF=indf_list_element_proc
        
        #replace PLACEHOLDERs 
        #with corresponding classes and properties
        
        #first, make the df with placeholder elements
        #and corresponding replacement classes
        #which we get from the entries with same ids
        #logic is that entityReference will have the same class 
        #as the corresponding name
        if(any(grepl("PLACEHOLDER$"
                    ,dF$class))){
            message("Placeholders found, fixing...")
            
            placeholder_df<-
                data_frame(rownum= grep("PLACEHOLDER$"
                                        ,dF$class)
                           ,id=dF$id[rownum]
                           ,to_replace=
                               dF$class[rownum])
            #get classes of elements w/o placeholders
            #which have ids of elements with placeholders
            placeholder_df$class<-
                dF[-placeholder_df$rownum,] %>%
                .[match(placeholder_df$id
                        ,.$id),] %>%
                .$class
            
            #case when match was not found -- avoid NA values
            placeholder_df$class[is.na(placeholder_df$class)]<-
                notfound
            
            ######netpath corner case of molecular interactions
            placeholder_df$class[grepl("Interaction"
                                       ,placeholder_df$id
                                       ,ignore.case =  TRUE)]<-
                "MolecularInteraction"
            ######netpath corner case of molecular interactions
                           
            #replace placeholder elements with actual classes
            dF$class<-
                dF$class %>%
                mgsub(pattern=placeholder_df$to_replace
                      ,replacement=placeholder_df$class
                      ,.
                      ,fixed=TRUE)

        }
        
        return(dF)
    }
######################################## internal_fix_placeholders ########################################
