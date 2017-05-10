######################################## internal_extract_last_path_element ########################################
internal_extract_last_path_element<-
    function(path_vector
             ,subdf
    ){
        #take out final path element -- 
        #name or db:id or positionStatus:sequencePosition
        #or (unlikely) some other element
        path_vector<-
            path_vector %>% 
            strsplit(split=":") %>% 
            unlist
       
        #find the property value for each element of the path vector 
        #(if there are more than one)
        value_vector<-
            sapply(path_vector
                   ,FUN=function(element){
                       #if more than one possible search terms
                       #split and search all
                       element<-
                           strsplit(element
                                    ,split="|"
                                    ,fixed=TRUE) %>%
                           unlist
                       temp_df<-
                           subdf %>%
                           #filter by property
                           .[.$property %in% element,] %>%
                           as.data.frame %>%
                           filter(property_value!="")
                       
                       if(nrow(temp_df)>0){
                           #mutate to get a line in the format
                           #property=value{class=class}
                           val<-
                               temp_df %>%
                               mutate(pr_prval=paste0(property
                                                      ,"="
                                                      ,property_value)
                                      ,pr_prval_cl=paste0(pr_prval
                                                          ,"{class="
                                                          ,class
                                                          ,"}")
                               ) %>%
                               #separate by ;
                               .$pr_prval_cl %>%
                               paste(collapse=";")
                           return(val)
                       } else {
                           return(NULL)
                       }
                       
                   }) %>%
            unlist %>%
            paste(collapse="::")
        if(value_vector==""){
            value_vector<-NA
        }
        
        return(value_vector)
    }
######################################## internal_extract_last_path_element ########################################
# ######################################## internal_extract_last_path_element ########################################
# internal_extract_last_path_element<-
#     function(path_vector
#              ,subdf
#     ){
#         #take out final path element -- name or db:id, respectively
#         #or (unlikely) some other element
#         
#         #if path vector implies search for db and id
#         #search for them in a pairwise manner
#         if (path_vector=="db:id"){
#             path_vector<-
#                 c("db"
#                   ,"id")
#         } else if (path_vector=="positionStatus:sequencePosition"){
#             path_vector<-
#                 c("positionStatus"
#                   ,"sequencePosition")
#         }
#         #find the property value for each element of the path vector 
#         #(if there are more than one)
#         value_vector<-
#             sapply(path_vector
#                    ,FUN=function(element){
#                        #if more than one possible search terms
#                        #split and search all
#                        element<-
#                            strsplit(element
#                                     ,split="|"
#                                     ,fixed=TRUE) %>%
#                            unlist
#                        temp_df<-
#                            subdf %>%
#                            #filter by property
#                            .[.$property %in% element,] %>%
#                            as.data.frame %>%
#                            filter(property_value!="")
#                        
#                        if(nrow(temp_df)>0){
#                            #mutate to get a line in the format
#                            #property=value{class=class}
#                            val<-
#                                temp_df %>%
#                                mutate(pr_prval=paste0(property
#                                                       ,"="
#                                                       ,property_value)
#                                       ,pr_prval_cl=paste0(pr_prval
#                                                           ,"{class="
#                                                           ,class
#                                                           ,"}")
#                                ) %>%
#                                #separate by ;
#                                .$pr_prval_cl %>%
#                                paste(collapse=";")
#                            return(val)
#                        } else {
#                            return(NULL)
#                        }
#                        
#                    }) %>%
#             unlist %>%
#             paste(collapse="::")
#         if(value_vector==""){
#             value_vector<-NA
#         }
# 
#         return(value_vector)
#     }
# ######################################## internal_extract_last_path_element ########################################
