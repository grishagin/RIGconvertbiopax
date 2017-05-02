######################################## extractBiopaxInteractions ########################################
extractBiopaxInteractions<-
    function(pw_id
             ,biopax
             ,biopax_source){
        
        #' @title
        #' Extract Pathway Interactions
        #' @description 
        #' Takes in a BioPAX object and returns a dataframe of all "interactions", one per row, for ONE specified pathway.
        #' @param pw_id Pathway id to process.
        #' @param biopax BioPAX object.
        #' @param biopax_source Name of the public source of the BioPAX object.
        
        #' @author 
        #' Ivan Grishagin
        
        message("Processing pathway with id ", pw_id)
        #get pathway names from pathway id/name dataframe
        pw_df<-
            listPathways(biopax)
        pw_name<-
            pw_df$name[pw_df$id==pw_id]
        #get pathway components dataframe
        #component ids
        pw_components_ids<-
            listPathwayComponents(biopax=biopax
                                  ,id=pw_id
                                  ,returnIDonly=TRUE
                                  ,biopaxlevel=biopax$biopaxlevel)
        
        if(is.null(pw_components_ids)){
            message("Pathway ",pw_id," does not have components, skipping.")
            return(NULL)
        }
        
        #get component instances INCLUDING all dependent entities
        dFrame<-
            selectInstances(biopax=biopax
                            ,id=pw_components_ids
                            ,includeReferencedInstances=TRUE)
        
        pathwayComponent_df<-
            adply(.data = pw_components_ids
                  ,.margins = 1
                  ,.fun = function(tempid){
                      temp_list<-
                          list(
                              pathwayComponent_id=tempid
                              ,pathwayComponent_class=
                                  filter(dFrame
                                         ,id==tempid) %>%
                                  .$class %>%
                                  unique
                              
                              ,controlType=
                                  filter(dFrame
                                         ,property=="controlType"
                                         ,id==tempid) %>%
                                  .$property_value
                              
                              ,catalysisDirection=
                                  filter(dFrame
                                         ,property=="catalysisDirection"
                                         ,id==tempid) %>%
                                  .$property_value
                              
                              ,spontaneous=
                                  filter(dFrame
                                         ,property=="spontaneous"
                                         ,id==tempid) %>%
                                  .$property_value
                              
                              ,controller_id=
                                  filter(dFrame
                                         ,property=="controller"
                                         ,id==tempid) %>%
                                  .$property_attr_value %>%
                                  striphash
                              
                              ,controller_class=NA
                              
                              ,controlled_id=
                                  filter(dFrame
                                         ,property=="controlled"
                                         ,id==tempid) %>%
                                  .$property_attr_value %>%
                                  striphash
                              
                              ,controlled_class=NA
                              ,noncontrol_id=NA
                              ,noncontrol_class=NA
                          ) 
                      #fill missing list elements with NA
                      temp_list[!unlist(lapply(temp_list,length))]<-NA
                      return(as.data.frame(temp_list))
                  }) %>%
            .[,-1]
       
        #add controller and controlled classes
        pathwayComponent_df$controller_class<-
            dFrame$class[match(pathwayComponent_df$controller_id
                               ,dFrame$id)]
        
        pathwayComponent_df$controlled_class<-
            dFrame$class[match(pathwayComponent_df$controlled_id
                               ,dFrame$id)]
        
        #add noncontrol components' ids and classes
        #which components are not controls?
        noncontrol_logi<-
            is.na(pathwayComponent_df$controller_id) &
            is.na(pathwayComponent_df$controlled_id)
        #add ids
        pathwayComponent_df$noncontrol_id[noncontrol_logi]<-
            pathwayComponent_df$pathwayComponent_id[noncontrol_logi]
        #add classes
        pathwayComponent_df$noncontrol_class[noncontrol_logi]<-
            pathwayComponent_df$pathwayComponent_class[noncontrol_logi]
                                          
        #get an evidence df at the pathway component level
        evidence_df<-
            adply(.data=pathwayComponent_df$pathwayComponent_id
                  ,.margins = 1
                  ,.fun = internal_find_entity_params
                  ,dFrame=dFrame
                  ,biopax_source=biopax_source
                  ,entity_type="pathwayComponent") %>%
            .[,-1]
        
        #get a df of controllers
        controller_df<-
            adply(.data=pathwayComponent_df$controller_id
                  ,.margins = 1
                  ,.fun = internal_find_entity_params
                  ,dFrame=dFrame
                  ,biopax_source=biopax_source
                  ,entity_type="controller") %>%
            .[,-1]
        
        #get a df of controlled entities
        controlled_df<-
            adply(.data=pathwayComponent_df$controlled_id
                  ,.margins = 1
                  ,.fun = internal_find_entity_params
                  ,dFrame=dFrame
                  ,biopax_source=biopax_source
                  ,entity_type="controlled") %>%
            .[,-1]
        
        #get a df of other_component (i.e. non-control-type) entities
        noncontrol_df<-
            adply(.data=pathwayComponent_df$noncontrol_id
                  ,.margins = 1
                  ,.fun = internal_find_entity_params
                  ,dFrame=dFrame
                  ,biopax_source=biopax_source
                  ,entity_type="noncontrol") %>%
            .[,-1]
        
        #merge dfs and write to file(s)
        result_df<-
            cbind.data.frame(pathway_id=pw_id
                             ,pathway_name=pw_name
                             ,source=biopax_source
                             ,pathwayComponent_df
                             ,evidence_df
                             ,controller_df
                             ,controlled_df
                             ,noncontrol_df)
        return(result_df)
    }
######################################## extractBiopaxInteractions ########################################

#' ######################################## extractBiopaxInteractions ########################################
#' extractBiopaxInteractions<-
#'     function(pw_id
#'              ,biopax
#'              ,biopax_source){
#'         
#'         #' @title
#'         #' Extract Pathway Interactions
#'         #' @description 
#'         #' Takes in a BioPAX object and returns a dataframe of all "interactions", one per row, for ONE specified pathway.
#'         #' @param pw_id Pathway id to process.
#'         #' @param biopax BioPAX object.
#'         #' @param biopax_source Name of the public source of the BioPAX object.
#'         
#'         #' @author 
#'         #' Ivan Grishagin
#'         
#'         message("Processing pathway with id ", pw_id)
#'         #get pathway names from pathway id/name dataframe
#'         pw_df<-
#'             listPathways(biopax)
#'         pw_name<-
#'             pw_df$name[pw_df$id==pw_id]
#'         #get pathway components dataframe
#'         #component ids
#'         pw_components_ids<-
#'             listPathwayComponents(biopax=biopax
#'                                   ,id=pw_id
#'                                   ,returnIDonly=TRUE
#'                                   ,biopaxlevel=biopax$biopaxlevel)
#'         
#'         if(is.null(pw_components_ids)){
#'             message("Pathway ",pw_id," does not have components, skipping.")
#'             return(NULL)
#'         }
#'         
#'         #get component instances INCLUDING all dependent entities
#'         dFrame<-
#'             selectInstances(biopax=biopax
#'                             ,id=pw_components_ids
#'                             ,includeReferencedInstances=TRUE)
#'         
#'         component_df<-
#'             adply(.data = pw_components_ids
#'                   ,.margins = 1
#'                   ,.fun = function(component_id){
#'                       temp_list<-
#'                           list(
#'                               component_id=component_id
#'                               ,component_class=
#'                                   filter(dFrame
#'                                          ,id==component_id) %>%
#'                                   .$class %>%
#'                                   unique
#'                               ,controlType=
#'                                   filter(dFrame
#'                                          ,property=="controlType"
#'                                          ,id==component_id) %>%
#'                                   .$property_value
#'                               
#'                               ,catalysisDirection=
#'                                   filter(dFrame
#'                                          ,property=="catalysisDirection"
#'                                          ,id==component_id) %>%
#'                                   .$property_value
#'                               
#'                               ,spontaneous=
#'                                   filter(dFrame
#'                                          ,property=="spontaneous"
#'                                          ,id==component_id) %>%
#'                                   .$property_value
#'                               
#'                               ,controller_id=
#'                                   filter(dFrame
#'                                          ,property=="controller"
#'                                          ,id==component_id) %>%
#'                                   .$property_attr_value %>%
#'                                   striphash
#'                               
#'                               ,controller_class=NA
#'                               
#'                               ,controlled_id=
#'                                   filter(dFrame
#'                                          ,property=="controlled"
#'                                          ,id==component_id) %>%
#'                                   .$property_attr_value %>%
#'                                   striphash
#'                               
#'                               ,controlled_class=NA
#'                               
#'                               ,non_ctrl_component=NA
#'                           ) 
#'                       #fill missing list elements with NA
#'                       temp_list[!unlist(lapply(temp_list,length))]<-NA
#'                       return(as.data.frame(temp_list))
#'                   }) %>%
#'             .[,-1]
#'         #non-control class components
#'         non_contr_comp<-
#'             apply(component_df[,c("controlType","controller_id","controlled_id")]
#'                   ,MARGIN = 1
#'                   ,FUN=function(row){
#'                       sum(is.na(row))
#'                   }) %>%
#'             as.logical
#'         component_df$non_ctrl_component<-
#'             non_contr_comp
#'         non_contr_comp[!non_contr_comp]<-
#'             NA
#'         #add controller and controlled classes
#'         component_df$controller_class<-
#'             dFrame$class[match(component_df$controller_id,dFrame$id)]
#'         
#'         component_df$controlled_class<-
#'             dFrame$class[match(component_df$controlled_id,dFrame$id)]
#'         
#'         #get a df of controllers
#'         controller_df<-
#'             adply(.data=component_df$controller_id
#'                   ,.margins = 1
#'                   ,.fun = internal_find_entity_params
#'                   ,dFrame=dFrame
#'                   ,biopax_source=biopax_source
#'                   ,entity_type="controller") %>%
#'             .[,-1]
#'         
#'         #get a df of controlled entities
#'         controlled_df<-
#'             adply(.data=component_df$controlled_id
#'                   ,.margins = 1
#'                   ,.fun = internal_find_entity_params
#'                   ,dFrame=dFrame
#'                   ,biopax_source=biopax_source
#'                   ,entity_type="controlled") %>%
#'             .[,-1]
#'         
#'         #get a df of other_component (i.e. non-control-type) entities
#'         non_ctrl_component_df<-
#'             adply(.data=component_df$component_id[non_contr_comp]
#'                   ,.margins = 1
#'                   ,.fun = internal_find_entity_params
#'                   ,dFrame=dFrame
#'                   ,biopax_source=biopax_source
#'                   ,entity_type="non_ctrl_component") %>%
#'             .[,-1]
#'         
#'         #merge dfs and write to file(s)
#'         result_df<-
#'             cbind.data.frame(pathway.ID=pw_id
#'                              ,pathway.Name=pw_name
#'                              ,Source=biopax_source
#'                              ,component_df
#'                              ,controller_df
#'                              ,controlled_df
#'                              ,non_ctrl_component_df)
#'         return(result_df)
#'     }
#' ######################################## extractBiopaxInteractions ########################################
