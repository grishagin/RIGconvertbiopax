######################################## internal_make_main_aux_df ########################################
internal_make_main_aux_df<-
    function(interactions_df){
        #makes main auxiliary dataframe based on pathway and component ids and classes
        #declare property attr values for different datatypes
        datatype_string<-
            "http://www.w3.org/2001/XMLSchema#string"
        datatype_boolean<-
            "http://www.w3.org/2001/XMLSchema#boolean"
        
        #declare key vectors
        pw_ids<-
            interactions_df$pathway_id
        pw_names<-
            interactions_df$pathway_name
        comp_ids<-
            interactions_df$pathwayComponent_id
        comp_classes<-
            interactions_df$pathwayComponent_class
        control_r_id<-
            interactions_df$controller_id
        control_d_id<-
            interactions_df$controlled_id
        control_type<-
            interactions_df$controlType
        catal_dir<-
            interactions_df$catalysisDirection
        spont<-
            interactions_df$spontaneous
        
        #declare list to store data tables
        dTable_list<-
            list()
        
        #declare and fill pathways biopax-style data table
        dTable_list$pw_ids<-
            data.frame(class="Pathway"
                       ,id=pw_ids
                       ,property="pathwayComponent"
                       ,property_attr="rdf:resource"
                       ,property_attr_value=comp_ids
                       ,property_value=""
            ) 
        #pathway names
        dTable_list$pw_names<-
            data.frame(class="Pathway"
                       ,id=pw_ids
                       ,property="name"
                       ,property_attr="rdf:datatype"
                       ,property_attr_value=datatype_string
                       ,property_value=pw_names
            ) 
        #declare and fill control class biopax-style data tables referring to controllers etc.
        #controller
        dTable_list$control_r_id<-
            data.frame(class=comp_classes
                       ,id=comp_ids
                       ,property="controller"
                       ,property_attr="rdf:resource"
                       ,property_attr_value=control_r_id
                       ,property_value=""
            )
        #controlled
        dTable_list$control_d_id<-
            data.frame(class=comp_classes
                       ,id=comp_ids
                       ,property="controlled"
                       ,property_attr="rdf:resource"
                       ,property_attr_value=control_d_id
                       ,property_value=""
            )
        #control types
        try(dTable_list$control_type<-
                data.frame(class=comp_classes
                           ,id=comp_ids
                           ,property="controlType"
                           ,property_attr="rdf:datatype"
                           ,property_attr_value=datatype_string
                           ,property_value=control_type
                ) 
            ,silent=TRUE)
        #catalysis direction
        try(dTable_list$catal_dir<-
                data.frame(class=comp_classes
                           ,id=comp_ids
                           ,property="catalysisDirection"
                           ,property_attr="rdf:datatype"
                           ,property_attr_value=datatype_string
                           ,property_value=catal_dir
                )
            ,silent=TRUE)
        #spontaneous?
        try(dTable_list$spont<-
                data.frame(class=comp_classes
                           ,id=comp_ids
                           ,property="spontaneous"
                           ,property_attr="rdf:datatype"
                           ,property_attr_value=datatype_boolean
                           ,property_value=spont
                )
            ,silent=TRUE)
        #bind, remove incomplete cases (i.e. with NA values)
        dTable<-
            dTable_list %>%
            do.call(rbind.data.frame
                    ,.) %>%
            .[complete.cases(.),] %>%
            unique
        rownames(dTable)<-NULL
        
        return(dTable)
    }
######################################## internal_make_main_aux_df ########################################

# ######################################## internal_make_main_aux_df ########################################
# internal_make_main_aux_df<-
#     function(interactions_df){
#         #makes main auxiliary dataframe based on pathway and component ids and classes
#         #declare property attr values for different datatypes
#         datatype_string<-
#             "http://www.w3.org/2001/XMLSchema#string"
#         datatype_boolean<-
#             "http://www.w3.org/2001/XMLSchema#boolean"
#         
#         #declare key vectors
#         pw_ids<-
#             interactions_df$pathway.ID
#         pw_names<-
#             interactions_df$pathway.Name
#         comp_ids<-
#             interactions_df$component_id
#         comp_classes<-
#             interactions_df$component_class
#         control_r_id<-
#             interactions_df$controller_id
#         control_d_id<-
#             interactions_df$controlled_id
#         control_type<-
#             interactions_df$controlType
#         catal_dir<-
#             interactions_df$catalysisDirection
#         spont<-
#             interactions_df$spontaneous
#         
#         #declare list to store data tables
#         dTable_list<-
#             list()
#         
#         #declare and fill pathways biopax-style data table
#         dTable_list$pw_ids<-
#             data.frame(class="Pathway"
#                        ,id=pw_ids
#                        ,property="pathwayComponent"
#                        ,property_attr="rdf:resource"
#                        ,property_attr_value=comp_ids
#                        ,property_value=""
#             ) 
#         #pathway names
#         dTable_list$pw_names<-
#             data.frame(class="Pathway"
#                        ,id=pw_ids
#                        ,property="name"
#                        ,property_attr="rdf:datatype"
#                        ,property_attr_value=datatype_string
#                        ,property_value=pw_names
#             ) 
#         #declare and fill control class biopax-style data tables referring to controllers etc.
#         #controller
#         dTable_list$control_r_id<-
#             data.frame(class=comp_classes
#                        ,id=comp_ids
#                        ,property="controller"
#                        ,property_attr="rdf:resource"
#                        ,property_attr_value=control_r_id
#                        ,property_value=""
#             )
#         #controlled
#         dTable_list$control_d_id<-
#             data.frame(class=comp_classes
#                        ,id=comp_ids
#                        ,property="controlled"
#                        ,property_attr="rdf:resource"
#                        ,property_attr_value=control_d_id
#                        ,property_value=""
#             )
#         #control types
#         try(dTable_list$control_type<-
#                 data.frame(class=comp_classes
#                            ,id=comp_ids
#                            ,property="controlType"
#                            ,property_attr="rdf:datatype"
#                            ,property_attr_value=datatype_string
#                            ,property_value=control_type
#                 ) 
#             ,silent=TRUE)
#         #catalysis direction
#         try(dTable_list$catal_dir<-
#                 data.frame(class=comp_classes
#                            ,id=comp_ids
#                            ,property="catalysisDirection"
#                            ,property_attr="rdf:datatype"
#                            ,property_attr_value=datatype_string
#                            ,property_value=catal_dir
#                 )
#             ,silent=TRUE)
#         #spontaneous?
#         try(dTable_list$spont<-
#                 data.frame(class=comp_classes
#                            ,id=comp_ids
#                            ,property="spontaneous"
#                            ,property_attr="rdf:datatype"
#                            ,property_attr_value=datatype_boolean
#                            ,property_value=spont
#                 )
#             ,silent=TRUE)
#         #bind, remove incomplete cases (i.e. with NA values)
#         dTable<-
#             dTable_list %>%
#             do.call(rbind.data.frame
#                     ,.) %>%
#             .[complete.cases(.),] %>%
#             unique
#         rownames(dTable)<-NULL
#         
#         return(dTable)
#     }
# ######################################## internal_make_main_aux_df ########################################
