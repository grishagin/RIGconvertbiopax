######################################## internal_find_entity_params ########################################
internal_find_entity_params<-
    function(entity_pav
             ,dFrame
             ,biopax_source
             ,entity_type=c("pathwayComponent","controller","controlled","noncontrol")
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
              ,"3participant"
              ,"4left"
              ,"5right"
            )
        path_keys[[2]]<-
            c("1NA"
              ,"2component"
            )
        
        path_keys[[3]]<-
            c("1name"
              ,"2entityReference_xref_dbid"
              ,"3xref_dbid"
              ,"4feature_modificationType_name"
              ,"4feature_modificationType_xref_dbid"
              ,"5feature_featureLocation_statpos"
              ,"5feature_sequenceIntervalBegin_statpos"
              ,"5feature_sequenceIntervalEnd_statpos"
              ,"6cellularLocation_name"
              ,"6cellularLocation_xref_dbid"
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
        #remove some properties based on entity_type
        #e.g. a biochemical reaction (left/right components) can't be a controller
        if(entity_type=="controller"){
            path_keys_vect<-
                path_keys_vect %>% 
                .[!grepl("product|participant|left|right"
                         ,.)]
        }else if(entity_type=="controlled"){
            path_keys_vect<-
                path_keys_vect %>% 
                .[!grepl("participant"
                         ,.)] %>% 
                c("evidence_evidenceCode_name"
                  ,"evidence_evidenceCode_xref_dbid"
                  ,.) 
        }else if(entity_type=="pathwayComponent"){
            #for high-level control component  
            #just take the evidence
            path_keys_vect<-
                c("evidence_evidenceCode_name"
                  ,"evidence_evidenceCode_xref_dbid")
        }
        ########################### prepare path keys vector
        
        ########################### actual path property elements
        path_elements_list<-
            list(
                product = "product"
                ,left = "left"
                ,right = "right"
                ,participant = "participant"
                ,component = "component"
                ,name = name_name
                ,entityReference = "entityReference"
                ,xref = c("xref"
                          ,"db:id")
                ,featureLocation = c("featureLocation"
                                     ,"positionStatus:sequencePosition")
                ,sequenceIntervalBegin = c("featureLocation"
                                           ,"sequenceIntervalBegin"
                                           ,"positionStatus:sequencePosition")
                ,sequenceIntervalEnd = c("featureLocation"
                                         ,"sequenceIntervalEnd"
                                         ,"positionStatus:sequencePosition")
                ,evidence = "evidence"
                ,evidenceCode = "evidenceCode"
                ,feature = "feature"
                ,modificationType = "modificationType"
                ,cellularLocation = "cellularLocation"
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

# ######################################## internal_find_entity_params ########################################
# internal_find_entity_params<-
#     function(entity_pav
#              ,dFrame
#              ,biopax_source
#              ,entity_type=c("controller","controlled","non_ctrl_component")
#     ){
#         if(nrow(dFrame)<1){
#             dFrame<-
#                 dFrame[1,]
#         }
#         #get entity class
#         entity_class<-
#             unique(dFrame$class[dFrame$id==entity_pav])
#         
#         if (length(entity_class)<1){
#             entity_class<-NA
#         }
#         
#         if (length(entity_class)!=1){
#             stop("There is a problem with entity class for entity "
#                  ,entity_pav
#                  ,". Found "
#                  ,length(entity_class)
#                  ," classes. Stopping.")
#         }
#         name_name<-"name|displayName|standardName|term"
#        
#         ########################### prepare path keys vector
#         path_keys<-list()
#         #create path keys -- to assemble desired path sequences
#         #numbers are added ONLY to ensure desired sequence of keys
#         #when assembled in one vector
#         path_keys[[1]]<-
#             c("1NA"
#               ,"2product"
#               ,"3left"
#               ,"4right"
#             )
#         path_keys[[2]]<-
#             c("1NA"
#               #,"memberPhysicalEntity" -- deprecated
#               ,"2complexcomponent"
#             )
#         
#         path_keys[[3]]<-
#             c("1name"
#               ,"2entityReference_xref_dbid"
#               ,"3xref_dbid"
#             )
#         
#         path_keys_vect<-
#             expand.grid(path_keys) %>%
#             apply(MARGIN=1
#                   ,paste0
#                   ,collapse="_") %>%
#             #remove NA values - these are to access
#             #1st level elements
#             gsub("NA_"
#                  ,""
#                  ,.) %>%
#             sort() %>%
#             #remove auxiliary digits
#             gsub("[[:digit:]]{1}"
#                  ,""
#                  ,.) 
#         ########################### prepare path keys vector
#         
#         ########################### actual path property elements
#         path_elements_list<-
#             list(
#                 product = "product"
#                 ,left = "left"
#                 ,right = "right"
#                 ,memberPhysicalEntity = "memberPhysicalEntity"
#                 ,complexcomponent = "component"
#                 ,name = name_name
#                 ,entityReference = "entityReference"
#                 ,xref = c("xref"
#                           ,"db:id")
#             )
#         
#         paths_list<-
#             path_keys_vect %>%
#             strsplit(split="_") %>%
#             lapply(FUN=function(element){
#                 path_elements_list[element] %>%
#                     unlist %>%
#                     stripNames
#             })
#         names(paths_list)<-
#             path_keys_vect
#         
#         ########################### actual path property elements
#         
#         
#         var_list<-
#             lapply(paths_list
#                    ,FUN=function(path,dFrame,entity_pav){
#                        internal_follow_path_extract_value(dFrame = dFrame
#                                                           ,pid = entity_pav
#                                                           ,path_vector=path)
#                    }
#                    ,dFrame=dFrame
#                    ,entity_pav=entity_pav)
#         
#         entity_df<-
#             as.data.frame(var_list)
#         colnames(entity_df)<-
#             paste(entity_type
#                   ,colnames(entity_df)
#                   ,sep="_")
#         
#         return(entity_df)
#     }
# ######################################## internal_find_entity_params ########################################
