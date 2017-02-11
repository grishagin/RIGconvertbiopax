######################################## pathway2RegulatoryGraph_Rancho ########################################
pathway2RegulatoryGraph_Rancho<-
    function (biopax
              ,pwid
              ,expandSubpathways = TRUE
              ,splitComplexMolecules = TRUE
              ,useIDasNodenames = FALSE
              ,verbose = FALSE
              ,returnGraph=TRUE){
        #' @title
        #' Pathway to Regulatory Graph
        #' @description 
        #' A different interpretation of a \code{pathway2RegulatoryGraph} function from rBiopaxParser package. 
        #' @details 
        #' Returns either a graph or a dataframe with control-type components.
        #' @param biopax BioPAX object.
        #' @param pwid Pathway ID.
        #' @param expandSubpathways Logical. Expand sub-pathways? Defaults to \code{TRUE}.
        #' @param splitComplexMolecules Logical. Split complex molecules into individual components? Defaults to \code{TRUE}.
        #' @param useIDasNodenames Logical. Use node ids as names? Defaults to \code{FALSE}.
        #' @param verbose Logical. Show all pertaining progress? Defaults to \code{FALSE}.
        #' @param returnGraph Logical. Return graph (as opposed to a dataframe of control-type components)? Defaults to \code{TRUE}.
        
        #' @author 
        #' Ivan Grishagin
        
        if (is.null(biopax) || 
            !("biopax" %in% class(biopax))) 
            stop("Error: pathway2RegulatoryGraph_Rancho: parameter biopax has to be of class biopax.")
        biopaxlevel<-
            biopax$biopaxlevel
        if (is.null(pwid) || 
            !("character" %in% class(pwid))) 
            stop("Error: pathway2RegulatoryGraph_Rancho: parameter pwid has to be of class character.")
        if (!require(graph)) {
            message(paste("This functions needs the graph library installed, albeit it cannot be found. Check out the installation instructions!" 
                          ,"\n"))
            return(NULL)
        }
        pwid<-
            unique(striphash(pwid))
        mygraph<-
            new(getClassDef("graphNEL"
                            ,package = "graph"))
        graph::edgemode(mygraph)<-"directed"
        pw_component_list<-
            listPathwayComponents(biopax
                                  ,pwid
                                  ,returnIDonly = T)
        if (length(pw_component_list) == 0) {
            warning("Pathway seems to have no pathway components")
            return(NULL)
        }
        pw_component_list<-
            selectInstances(biopax
                            ,id = pw_component_list
                            ,includeReferencedInstances = TRUE
                            ,returnCopy = TRUE)
        pw_component_list$property = tolower(pw_component_list$property)
        setkeyv(pw_component_list, cols = c("id"
                                            ,"class"
                                            ,"property"))
        pw_controls<-
            pw_component_list[tolower(pw_component_list$class) %chin% 
                                  c("control"
                                    ,"catalysis"
                                    ,"modulation"
                                    ,"templatereactionregulation"), ]
        if (length(pw_controls$id) == 0) {
            warning("warning: pathway2RegulatoryGraph_Rancho: supplied graph has no regulatory pathway components. Returning NULL.")
            return(NULL)
        } else {
            if (verbose) {
                message(paste("Found"
                              ,length(unique(pw_controls$id))
                              ,"control-type pathway components. Putting them together..."))
            }
        }
        #dataframe to store all controllers and controlleds
        control_df<-NULL
        for (i in sort(unique(pw_controls$id))) {
            instance<-
                pw_controls[id == i, ]
            if (biopax$biopaxlevel == 2) {
                type<-
                    as.character(instance[property == "control-type"]$property_value)
            }
            if (biopax$biopaxlevel == 3) {
                type<-
                    as.character(instance[property == "controltype"]$property_value)
            }
            if (length(type) == 0) {
                if (tolower(as.character(instance[1, class])) == "catalysis") {
                    type<-
                        "ACTIVATION"
                }
                else {
                    next
                }
            }
            if (!grepl("activation"
                       ,type
                       ,ignore.case = TRUE) & 
                !grepl("inhibition"
                       ,type
                       ,ignore.case = TRUE)) {
                next
            }
            controller_ids<-
                striphash(as.character(unique(instance[property == "controller"]$property_attr_value)))
            controllers<-
                NA
            for (i2 in controller_ids) {
                c_instance<-
                    pw_component_list[id == i2, ]
                if (biopax$biopaxlevel == 2) {
                    c_instance<-
                        pw_component_list[id == striphash(c_instance[property == "physical-entity"]$property_attr_value), ]
                }
                if (splitComplexMolecules & 
                    any(isOfClass(c_instance 
                                  ,"complex"))) {
                    if (useIDasNodenames) {
                        controllers<-
                            c(controllers
                              ,as.character(splitComplex(pw_component_list 
                                                         ,i2
                                                         ,returnIDonly = T
                                                         ,biopaxlevel = biopaxlevel)))
                    } else {
                        controllers<-
                            c(controllers
                              ,as.character(splitComplex(pw_component_list
                                                         ,i2
                                                         ,biopaxlevel = biopaxlevel)$name))
                    }
                } else {
                    if (useIDasNodenames) {
                        controllers<-
                            c(controllers
                              ,c_instance$id[1])
                    } else {
                        controllers<-
                            c(controllers
                              ,getInstanceProperty(pw_component_list
                                                   ,c_instance$id[1]
                                                   ,biopaxlevel = biopaxlevel))
                    }
                }
            }
            
            controlled_ids<-
                striphash(as.character(unique(instance[property == "controlled"]$property_attr_value)))
            controlleds<-
                NA
            for (i2 in controlled_ids) {
                c_instance<-
                    pw_component_list[id == i2, ]
                if (any(isOfClass(c_instance
                                  ,c("conversion")
                                  ,considerInheritance = TRUE)) || 
                    any(isOfClass(c_instance, c("templatereaction")))) {
                    leftrights<-
                        striphash(c_instance[property == "left" | 
                                                 property == "right" | 
                                                 property == "product"]$property_attr_value)
                    for (i3 in leftrights) {
                        leftrights_instance <-
                            pw_component_list[id == i3, ]
                        if (biopax$biopaxlevel == 2) {
                            leftrights_instance<-
                                pw_component_list[id == striphash(leftrights_instance[property == "physical-entity"]$property_attr_value), ]
                        }
                        if (splitComplexMolecules & 
                            any(isOfClass(leftrights_instance 
                                          ,"complex"))) {
                            if (useIDasNodenames) {
                                controlleds<-
                                    c(controlleds, as.character(splitComplex(pw_component_list
                                                                             ,i3
                                                                             ,returnIDonly = T
                                                                             ,biopaxlevel = biopaxlevel)))
                            } else {
                                controlleds<-
                                    c(controlleds, as.character(splitComplex(pw_component_list
                                                                             ,i3
                                                                             ,biopaxlevel = biopaxlevel)$name))
                            }
                        } else {
                            if (useIDasNodenames) {
                                controlleds<-
                                    c(controlleds
                                      ,as.character(leftrights_instance[1]$id))
                            } else {
                                controlleds<-
                                    c(controlleds
                                      ,getInstanceProperty(pw_component_list
                                                           ,leftrights_instance[1]$id
                                                           ,biopaxlevel = biopaxlevel))
                            }
                        }
                    }
                }
            }
            controllers<-
                striphash(unique(controllers))
            controlleds<-
                striphash(unique(controlleds))
            controllers<-
                controllers[!is.na(controllers) & 
                                !is.null(controllers) & 
                                nchar(controllers) > 0]
            controlleds<-
                controlleds[!is.na(controlleds) & 
                                !is.null(controlleds) & 
                                nchar(controlleds) > 0]
            
            controllers<-
                sort(controllers)
            controlleds<-
                sort(controlleds)

            if (length(controllers) == 0 | length(controlleds) == 0) {
                next
            }
           
            if (tolower(type) == "activation") {
                weights = 1
            } else {
                weights = -1
            }
            #add all controllers and controlleds combinations to df
            control_df<-
                expand.grid(controllers=controllers
                            ,controlleds=controlleds
                            ,stringsAsFactors=FALSE) %>%
                cbind.data.frame(weights=weights) %>%
                rbind.data.frame(control_df)
            
        }
        if(!is.null(control_df)){
            return(control_df)
        }
        control_df<-
            control_df %>%
            unique %>%
            filter(controllers!=controlleds) %>%
            arrange(controllers
                    ,controlleds)
        
        if(returnGraph){
            newnodes<-
                unique(c(control_df$controllers
                         ,control_df$controlleds))
            
            mygraph<- 
                graph::addNode(newnodes
                               ,mygraph)
            for(rindex in 1:nrow(control_df)){
                mygraph<-
                    graph::addEdge(control_df$controllers[rindex]
                                   ,control_df$controlleds[rindex]
                                   ,mygraph
                                   ,weights = weights)
            }
            
            return(mygraph)
        } else {
            return(control_df)
        }
        
       
    }
######################################## pathway2RegulatoryGraph_Rancho ########################################
