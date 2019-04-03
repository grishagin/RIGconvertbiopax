######################################## internal_string_to_df_inner ########################################
internal_string_to_df_inner<-
    function(string
             ,el_id=""
             ,prop_vect){
        #first, split by inner separator, which is a pipe
        #it's for multiple same-level xref entries
        stringvect<-
            string %>%
            strsplit(split="|"
                     ,fixed=TRUE) %>%
            unlist
        #fabricate index to distinguish
        #different entities 
        if(length(stringvect)>1){
            #if more than one pipe found,
            #prepare property attr values in advance
            #so basically we'd have the same el_id
            #refer to xref values with diff indices
            pavs<-
                paste(el_id
                      ,1:length(stringvect)
                      ,sep = "-")
        } else {
            pavs<-
                el_id
        }
        
        result_df<-
            stringvect %>%
            seq_along %>%
            lapply(function(lindex){
                temp_df<-
                    internal_string_to_df_inner_inner(string=stringvect[[lindex]]
                                                      ,el_id=el_id
                                                      ,pav=pavs[[lindex]]
                                                      ,prop_vect=prop_vect)
                
                return(temp_df)
            }) %>%
            do.call(rbind.data.frame
                    ,.)
        return(result_df)
    }
######################################## internal_string_to_df_inner ########################################
