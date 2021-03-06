######################################## internal_follow_path_extract_value ########################################
internal_follow_path_extract_value<-
    function(dFrame
             ,pid
             ,path_vector
             ,cpxlvl=1
             ,lvl=1
    ){
        #extract a df by id
        subdf<-
            dFrame[dFrame$id==pid,]
        
        #get a vector of referenced ids
        pav_vector<-
            subdf %>%
            .[.$property==path_vector[1],] %>%
            .$property_attr_value %>%
            striphash
        
        #perform a check for a complex with a name, but no components
        #no pavs found
        if(length(pav_vector)<1 &
           #current property that is being looked for is "component"
           path_vector[1]=="component" &
           #and there's a complex amongst the classes of the corresponding df
           "Complex" %in% subdf$class){
            #in that case, move to the next path property (i.e. name or xref)
            path_vector<-
                path_vector[-1]
        }
        
        #in case it is the end of the path vector, convert to string with name and class
        if(length(path_vector)==1){
            to_return<-
                internal_extract_last_path_element(path_vector=path_vector
                                                   ,subdf=subdf)
            #from if-clause return value vector
            return(to_return)
        } 
        
        #if those ids cannot be found - abort
        if (sum(pav_vector %in% dFrame$id)==0){
            #message("internal_follow_path_extract_value: pav_vector not found in ids.")
            return(NA)
        }
        #go through the vector of referenced ids recursively
        #add different separators between different levels
        #if such level is "xref", then just use one pipe to 
        #not separate them apart too much
        
        #if(path_vector[1]=="xref"){
        if(path_vector[1] %in% internal_banned_props()){
            sep<-
                "|"
        }else if(path_vector[1]=="component" & 
                 cpxlvl>1){
            sep<-
                paste0("|",cpxlvl,"|")
        }else{
            #assign each level its separator ||, |-|, etc.
            sep<-
                paste0("|",rep("-",lvl-1),"|")
        }
        
        #find present component's name or xref 
        #(depending on what the next path element is)
        #basically, if current path element is component 
        #and current dFrame contains component properties
        #there should also be a complex name or xref
        #which we can extract
        #like a side branch from  program
        if(path_vector[1]=="component" &
           "component" %in% dFrame$property &
           #only for NOT top-level component
           #that one has its name in the preceding *blah-blah-_name* column
           #I know, it's confusing
           cpxlvl>1){
            complex_name<-
                internal_follow_path_extract_value(dFrame=dFrame
                                                   ,pid=pid
                                                   #basically, move onto the next property
                                                   ,path_vector=path_vector[-1]
                                                   ,cpxlvl=1) %>%
                #complex name is in parentheses
                paste0("(",.,")")
        }
        
        to_return<-
            sapply(pav_vector
                   ,FUN=function(pav){
                       #in general case, take first item off of the path vector
                       path_vector_new<-
                           path_vector[-1]
                       #check if there's a complex within a complex
                       #if the current path vector is component...
                       if(path_vector[1]=="component"){
                           #get a df for corresponding pav
                           tempdf<-
                               dFrame[dFrame$id==pav,]
                           
                           #check class of the returned value(s)
                           #if it's complex, it means there's another complex
                           #coming up, i.e. complex in complex, so need to repeat
                           #"component"path property
                           #and change separator
                           
                           if("complex" %in% tolower(tempdf$class)){
                               #keep looking for components
                               path_vector_new<-
                                   path_vector
                               #change the level
                               cpxlvl<-
                                   cpxlvl+1
                           }
                           
                       }
                       val<-
                           internal_follow_path_extract_value(dFrame=dFrame
                                                              ,pid=pav
                                                              ,path_vector=path_vector_new
                                                              ,cpxlvl=cpxlvl
                                                              ,lvl=lvl+1
                           )
                       return(val)
                   }) %>%
            unlist %>%
            paste(collapse=sep) 
        
        #if a complex was just processed,
        #wrap the discovered values in a special pair of parentheses
        #they are required later, to be able to "unwrap" the complex
        if(exists("complex_name"
                  ,envir=environment()
                  ,inherits=FALSE) & 
           !is.na(to_return) & 
           !is.null(to_return)){
            to_return<-
                to_return %>%
                paste0("cplx"
                       ,complex_name
                       ,"<<\uff5f"
                       ,.
                       ,"\uff60>>cplx"
                )
        }
        return(to_return)
        
    }
######################################## internal_follow_path_extract_value ########################################
# ######################################## internal_follow_path_extract_value ########################################
# internal_follow_path_extract_value<-
#     function(dFrame
#              ,pid
#              ,path_vector
#              ,cpxlvl=1
#              ,lvl=1
#              #,drill_down=FALSE
#     ){
#         #extract a df by id
#         subdf<-
#             dFrame[dFrame$id==pid,]
#         
#         #get a vector of referenced ids
#         pav_vector<-
#             subdf %>%
#             .[.$property==path_vector[1],] %>%
#             .$property_attr_value %>%
#             striphash
#         
#         #perform a check for a complex with a name, but no components
#         if(length(pav_vector)<1 &
#            path_vector[1]=="component"){
#             #in that case, move to the next path property (i.e. name or xref)
#             path_vector<-
#                 path_vector[-1]
#         }
#         
#         #in case it is the end of the path vector, convert to string with name and class
#         if(length(path_vector)==1){
#             to_return<-
#                 internal_extract_last_path_element(path_vector=path_vector
#                                                    ,subdf=subdf)
#             #from if-clause return value vector
#             return(to_return)
#         } 
#         
#         #deprecated
#         #if the value has not been found, repeat the search with 
#         #the next property
#         #and then the next one - i.e skip levels, in case in some biopax files
#         #some levels can be missing
#         #by default - do not, i.e. the first (top) level HAS to be present
#         # if(drill_down){
#         #     while(length(pav_vector)==0 & 
#         #           length(path_vector>0)){
#         #         path_vector<-path_vector[-1]
#         #         pav_vector<-
#         #             subdf %>%
#         #             .[.$property==path_vector[1],] %>%
#         #             .$property_attr_value %>%
#         #             striphash
#         #     }
#         # }
#         
#         #if those ids cannot be found - abort
#         if (sum(pav_vector %in% dFrame$id)==0){
#             #message("internal_follow_path_extract_value: pav_vector not found in ids.")
#             return(NA)
#         }
#         #go through the vector of referenced ids recursively
#         #add different separators between different levels
#         #if such level is "xref", then just use one pipe to 
#         #not separate them apart too much
#         
#         if(path_vector[1]=="xref"){
#             sep<-
#                 "|"
#         }else if(path_vector[1]=="component" & 
#                  cpxlvl>1){
#             sep<-
#                 paste0("|",cpxlvl,"|")
#         }else{
#             #assign each level its separator ||, |-|, etc.
#             sep<-
#                 paste0("|",rep("-",lvl-1),"|")
#         }
#         
#         #find present component's name or xref 
#         #(depending on what the next path element is)
#         #basically, if current path element is component 
#         #and current dFrame contains component properties
#         #there should also be a complex name or xref
#         #which we can extract
#         #like a side branch from  program
#         if(path_vector[1]=="component" &
#            "component" %in% dFrame$property &
#            #only for NOT top-level component
#            #that one has its name in the corresponding Name column
#            #I know, it's confusing
#            cpxlvl>1){
#             complex_name<-
#                 internal_follow_path_extract_value(dFrame=dFrame
#                                                    ,pid=pid
#                                                    #basically, path_vector[2] is either a name or xref
#                                                    ,path_vector=path_vector[2]
#                                                    ,cpxlvl=1
#                                                    #,drill_down=FALSE
#                 ) %>%
#                 #complex name is in vertical "brackets"
#                 paste0("(",.,")")
#         }
#         
#         to_return<-
#             sapply(pav_vector
#                    ,FUN=function(pav){
#                        #in general case, take first item off of the path vector
#                        path_vector_new<-
#                            path_vector[-1]
#                        #check if there's a complex within a complex
#                        #if the current path vector is component...
#                        if(path_vector[1]=="component"){
#                            #get a df for corresponding pav
#                            tempdf<-
#                                dFrame[dFrame$id==pav,]
#                            
#                            #check class of the returned value(s)
#                            #if it's complex, it means there's another complex
#                            #coming up, i.e. complex in complex, so need to repeat
#                            #"component"path property
#                            #and change separator
#                            
#                            if("complex" %in% tolower(tempdf$class)){
#                                #keep looking for components
#                                path_vector_new<-
#                                    path_vector
#                                #change the level
#                                cpxlvl<-
#                                    cpxlvl+1
#                            }
#                            
#                        }
#                        val<-
#                            internal_follow_path_extract_value(dFrame=dFrame
#                                                               ,pid=pav
#                                                               ,path_vector=path_vector_new
#                                                               ,cpxlvl=cpxlvl
#                                                               ,lvl=lvl+1
#                                                               #,drill_down=FALSE
#                            )
#                        return(val)
#                    }) %>%
#             unlist %>%
#             paste(collapse=sep) 
#         
#         if(exists("complex_name"
#                   ,envir=environment()
#                   ,inherits=FALSE) & 
#            !is.na(to_return) & 
#            !is.null(to_return)){
#             to_return<-
#                 to_return %>%
#                 paste0("cplx"
#                        ,complex_name
#                        ,"<<\uff5f"
#                        ,.
#                        ,"\uff60>>cplx"
#                 )
#         }
#         return(to_return)
#         
#     }
# ######################################## internal_follow_path_extract_value ########################################
