######################################## internal_cplx_string_to_df ########################################
internal_cplx_string_to_df<-
    function(string
             ,el_id
             ,prop_vect){
        #declare the split between complex components
        split<-
            "\\|[[:digit:]]\\|"
        
        #build a df from complex components
        cplx_df<-
            internal_cplx_split(string) %>%
            split_cols_lengthen_df(colsToSplit = "comp"
                                   ,patternToSplit = split) %>%
            as.data.table %>% 
            #remove all components that are NA
            .[!is.na(comp)] %>% 
            #remove all components that refer to complex ids that are not in id column
            .[!(grepl("cplx\\."
                      ,comp) &
                    !comp %in% cplx.id)]
        if(nrow(cplx_df)==0){
            return(NULL)
        }
        
        #add all other columns required for biopax table
        cplx_df$class<-
            "Complex"
        cplx_df$id<-
            cplx_df$cplx.id
        cplx_df$property<-
            "component"
        cplx_df$property_attr<-
            "rdf:resource"
        cplx_df$property_attr_value<-
            ""
        cplx_df$property_value<-
            ""
        
        #comp column has either ids of complexes inside complexes
        #or names of entities
        #hence for those values, that are complex ids --
        #-- add them to property attr (reference) column
        cplx_df[comp %in% cplx.id]$property_attr_value<-
            cplx_df[comp %in% cplx.id]$comp
        
        #and for those values, that are NOT complex ids -- 
        #make up property attr values by addin _cc to every id
        #as well as a number from 1 to number of referenced components
        cplx_df[!comp %in% cplx.id]$property_attr_value<-
            cplx_df[!comp %in% cplx.id
                    ,.(property_attr_value=paste0(cplx.id
                                                  ,"_cc"
                                                  ,1:length(property_attr_value)))
                    ,by=cplx.id]$property_attr_value 
        
        #transform the protein and complex names from strings into
        #dataframes suitable for biopax data table
        split_names_df<-
            lapply(1:nrow(cplx_df)
                   ,FUN=function(rndex){
                       rbind.data.frame(internal_string_to_df_inner(string = cplx_df$cplx.name[rndex]
                                                                    ,el_id = cplx_df$cplx.id[rndex]
                                                                    ,prop_vect = prop_vect)
                                        ,internal_string_to_df_inner(string = cplx_df$comp[rndex]
                                                                     ,el_id = cplx_df$property_attr_value[rndex]
                                                                     ,prop_vect = prop_vect)
                       )
                   }) %>%
            do.call(rbind.data.frame
                    ,.) %>%
            as.data.table %>% 
            #remove duplicate rows
            unique
        
        
        result_df<-
            #select just the necessary columns
            cplx_df[,.(class
                       ,id
                       ,property 
                       ,property_attr                     
                       ,property_attr_value
                       ,property_value)][property_value!="toremove"] %>% 
            #bind with "names" dataframe
            rbind(split_names_df)
        #modify ids and attr_values  by adding general element ids
        result_df$id<-
            paste(el_id
                  ,result_df$id
                  ,sep="-")
        result_df[property_attr=="rdf:resource"]$property_attr_value<-
            paste(el_id
                  ,result_df[property_attr=="rdf:resource"]$property_attr_value
                  ,sep="-")
        result_df[grep("PLACEHOLDER"
                       ,class)]$class<-
            paste(el_id
                  ,result_df[grep("PLACEHOLDER"
                                  ,class)]$class
                  ,sep="-")
        
        return(result_df)
    }
######################################## internal_cplx_string_to_df ########################################
