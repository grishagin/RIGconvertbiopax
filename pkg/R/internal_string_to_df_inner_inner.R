######################################## internal_string_to_df_inner_inner ########################################
internal_string_to_df_inner_inner<-
    function(string
             ,el_id=""
             ,pav=""
             ,prop_vect){
        require(stringr)
        pattern = "([[:alpha:]]*?)=(.*?)\\{class=(.*?)\\}"
        #search for a given pattern in the string
        
        dFrame<-
            string %>%
            stringr::str_match_all(pattern) %>%
            .[[1]] %>%
            as.data.frame
        
        #if resultant dFrame does not have rows, return null
        if(nrow(dFrame)<1){
            return(NULL)
        }
        
        #... and take only pertaining columns of the hits
        #these numbers are hard-coded due to the pattern is known 
        #and expected to be one and only
        dFrame<-dFrame[,c(2:4)]
        
        #again, pattern is known and that's what it is
        colnames(dFrame)<-
            c("property"
              ,"property_value"
              ,"class")
        
        #exclude all NA values
        dFrame[dFrame=="NA"]<-NA
        dFrame<-
            dFrame %>% 
            filter(!is.na(property_value))
        
        #if resultant dFrame does not have rows, return null
        if(nrow(dFrame)<1){
            return(NULL)
        }
        
        #if it's a string with xrefs -- add corresponding 
        #joining lines (see below) and change ids accordingly
        #to account for appropriate id nesting
        if("entityReference" %in% prop_vect){
            #since the class to which the xref belongs is not known
            #use a PLACEHOLDER instead
            aux_df<-
                data.frame(class=c(paste0(el_id
                                          ,"PLACEHOLDER")
                                   ,paste0(el_id
                                           ,"PLACEHOLDERReference"))
                           ,id=c(el_id
                                 ,paste0(pav
                                         ,"_x"))
                           ,property=c("entityReference"
                                       ,"xref")
                           ,property_attr="rdf:resource"
                           ,property_attr_value=c(paste0(pav
                                                         ,"_x")
                                                  ,paste0(pav
                                                          ,"_xx"))
                           ,property_value="")
            
            el_id<-
                paste0(pav
                       ,"_xx")
            
        } else if ("xref" %in% prop_vect) {
            aux_df<-
                data.frame(class=c(paste0(el_id
                                          ,"PLACEHOLDER"))
                ,id=c(el_id)
                ,property=c("xref")
                ,property_attr="rdf:resource"
                ,property_attr_value=c(paste0(pav
                                              ,"_xr"))
                ,property_value="")
            
            el_id<-
                paste0(pav
                       ,"_xr")
        } else {
            #in case, there's no xref references,
            #there's no need for auxiliary df
            aux_df<-NULL
        }

        #build element df
        dFrame$id<-el_id 
        dFrame$property_attr<-"rdf:datatype"
        dFrame$property_attr_value<-"http://www.w3.org/2001/XMLSchema#string"
        
        #change order
        dFrame<-
            dFrame %>%
            dplyr::select(class
                          ,id
                          ,property
                          ,property_attr
                          ,property_attr_value
                          ,property_value) %>%
            #add that auxiliary dataframe
            rbind.data.frame(aux_df)
        
        return(dFrame)
    }
######################################## internal_string_to_df_inner_inner ########################################
