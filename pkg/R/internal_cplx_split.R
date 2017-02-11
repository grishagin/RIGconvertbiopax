######################################## internal_cplx_split ########################################
internal_cplx_split<-
    function(string
             ,lvl=1){
        require(stringr)
        require(qdap)
        #define lowest level complex pattern
        #i.e. cplx(NAME)<<[COMPONENTS THAT ARE NOT COMPLEX]>>cplx
        pattern<-
            #"cplx\\(([^\\)\\)]*)\\)<<\\[([^<<\\[\\]>>]*)\\]>>cplx"
            #"cplx¦([^¦]*)?¦<<\\[([^<<\\[\\]>>]*)\\]>>cplx"
            #"cplx\\((.*?)\\)<<\\[((?!cplx).*)\\]>>cplx"
            "cplx\\(([^\uff5f]*?)\\)<<\uff5f([^\uff5f\uff60]*)\uff60>>cplx"
        
        #dataframe with all complexes in one complex
        cplx_df<-
            string %>%
            str_match_all(pattern) %>%
            #stri_match_all(regex=pattern) %>%
            .[[1]] %>%
            as.data.frame 
        if(nrow(cplx_df)<1){
            message("internal_cplx_split: cplx_df has 0 rows for string: \n"
                    ,string)
        }
        
        #add discovered lowest level complex ids
        cplx_df$cplx.id<-
            paste("cplx",lvl,1:nrow(cplx_df)
                  ,sep=".")
        
        #replace those complexes in original string
        #with their ids
        string<-
            string %>%
            mgsub(cplx_df$V1
                  ,cplx_df$cplx.id
                  ,.)
        #take pertaining columns
        cplx_df<-
            cplx_df %>%
            dplyr::select(cplx.id
                          ,cplx.name=V2
                          ,comp=V3)
        
        #if there're more complexes in the new string
        #repeat the procedure recursively, append results
        if(length(grep("\uff60>>cplx"
                       ,string))>0){
            cplx_new_df<-
                internal_cplx_split(string
                           ,lvl=lvl+1)
            cplx_df<-
                rbind.data.frame(cplx_new_df
                                 ,cplx_df)
        }
        
        return(cplx_df)
    }
######################################## internal_cplx_split ########################################
