######################################## internal_find_entity_params ########################################
internal_find_entity_params<-
    function(entity_pav
             ,dFrame
             ,biopax_source
             ,entity_type=c("controller","controlled","non_ctrl_component")
    ){
        if(nrow(dFrame)<1){
            dFrame<-
                dFrame[1,]
        }
        #get entity class
        entity_class<-
            unique(dFrame$class[dFrame$id==entity_pav])
        
        if (length(entity_class)<1){
            entity_class<-NA
        }
        
        if (length(entity_class)!=1){
            stop("There is a problem with entity class for entity "
                 ,entity_pav
                 ,". Found "
                 ,length(entity_class)
                 ," classes. Stopping.")
        }
        name_name<-"name|displayName|standardName|term"
       
        ########################### prepare path keys vector
        path_keys<-list()
        #create path keys -- to assemble desired path sequences
        #numbers are added ONLY to ensure desired sequence of keys
        #when assembled in one vector
        path_keys[[1]]<-
            c("1NA"
              ,"2product"
              ,"3left"
              ,"4right"
            )
        path_keys[[2]]<-
            c("1NA"
              #,"memberPhysicalEntity" -- deprecated
              ,"2complexcomponent"
            )
        
        path_keys[[3]]<-
            c("1name"
              ,"2entityReference_xref_dbid"
              ,"3xref_dbid"
            )
        
        path_keys_vect<-
            expand.grid(path_keys) %>%
            apply(MARGIN=1
                  ,paste0
                  ,collapse="_") %>%
            #remove NA values - these are to access
            #1st level elements
            gsub("NA_"
                 ,""
                 ,.) %>%
            sort() %>%
            #remove auxiliary digits
            gsub("[[:digit:]]{1}"
                 ,""
                 ,.) 
        ########################### prepare path keys vector
        
        ########################### actual path property elements
        path_elements_list<-
            list(
                product = "product"
                ,left = "left"
                ,right = "right"
                ,memberPhysicalEntity = "memberPhysicalEntity"
                ,complexcomponent = "component"
                ,name = name_name
                ,entityReference = "entityReference"
                ,xref = c("xref"
                          ,"db:id")
            )
        
        paths_list<-
            path_keys_vect %>%
            strsplit(split="_") %>%
            lapply(FUN=function(element){
                path_elements_list[element] %>%
                    unlist %>%
                    stripNames
            })
        names(paths_list)<-
            path_keys_vect
        
        ########################### actual path property elements
        
        
        var_list<-
            lapply(paths_list
                   ,FUN=function(path,dFrame,entity_pav){
                       internal_follow_path_extract_value(dFrame = dFrame
                                                          ,pid = entity_pav
                                                          ,path_vector=path)
                   }
                   ,dFrame=dFrame
                   ,entity_pav=entity_pav)
        
        entity_df<-
            as.data.frame(var_list)
        colnames(entity_df)<-
            paste(entity_type
                  ,colnames(entity_df)
                  ,sep="_")
        
        return(entity_df)
    }
######################################## internal_find_entity_params ########################################
