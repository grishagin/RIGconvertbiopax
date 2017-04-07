######################################## internal_string_to_df ########################################
internal_string_to_df<-
    function(string
             ,el_id
             ,prop_vect
             ,split="||"){
        #function to split same level elements 
        #separated by split
        
        #split the vector of string elements by split pattern
        splitlist<-
            string %>%
            strsplit(split = split
                     ,fixed = TRUE) %>%
            unlist
        
        if(length(splitlist)<1){
            return(NULL)
        } else if(length(prop_vect[!prop_vect %in% RIGconvertbiopax:::internal_banned_props()])>0){
            
            #add names by combining the supplied element name
            #add list index
            names(splitlist)<-
                paste(el_id
                      ,1:length(splitlist)
                      ,sep="-"
                )
        } else {
            names(splitlist)<-
                el_id
        }
        
        #walk along the list of strings
        #and split these strings accordingly
        result<-
            splitlist %>%
            seq_along %>%
            #for each list element...
            lapply(X=.
                   ,splitlist=splitlist
                   ,prop_vect=prop_vect
                   ,FUN=function(lindex
                                 ,splitlist
                                 ,prop_vect){
                       #get list element...
                       el<-
                           splitlist[lindex]
                       #... and its id is equal to its name
                       el_id<-
                           names(splitlist)[lindex]
                       
                       #exclude NA and "NA"
                       if(is.na(el)){
                           return(NULL)
                       } else if(el=="NA"){
                           return(NULL)
                       }
                       
                       #try also splitting on another pattern
                       #if it is in the string
                       if(any(grepl("|-|"
                                    ,el
                                    ,fixed=TRUE))){
                           result<-
                               RIGconvertbiopax:::internal_string_to_df(string=el
                                                                        ,el_id=el_id
                                                                        ,prop_vect=prop_vect
                                                                        ,split="|-|")
                       } else if(length(grep("\uff60>>cplx"
                                             ,el))>0){
                           #if it's a complex, 
                           #launch complex processing function
                           result<-
                               RIGconvertbiopax:::internal_cplx_string_to_df(string=el
                                                                             ,el_id=el_id
                                                                             ,prop_vect=prop_vect)
                       } else {
                           #otherwise launch the general string processing function
                           result<-
                               RIGconvertbiopax:::internal_string_to_df_inner(string=el
                                                                              ,el_id=el_id
                                                                              ,prop_vect=prop_vect)
                       }
                       return(result)
                   }) %>%
            #bind them into a single dataframe
            do.call(rbind.data.frame
                    ,.) 
        return(result)
    }
######################################## internal_string_to_df ########################################
