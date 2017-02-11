######################################## MAIN_biopax_comparison ########################################
MAIN_biopax_comparison<-
    function(pwid_to_compare=c("all"
                               ,"none")
             ,pwid_to_plot=c("none"
                             ,"same"
                             ,"test5")
             ,new_biopax
             ,orig_biopax
             ,new_biopax_tag="_new"
             ,orig_biopax_tag="_orig"
             ,verbose = TRUE
             ,pw_df_path="default"
             ,source_name=NULL
             
             
    ){
        #' @title
        #' MAIN -- Compare BioPAX Objects
        #' @description 
        #' Compare pathways in two BioPAX objects. 
        #' @details 
        #' Allows to compare control components
        #' @param pwid_to_compare Which pathway IDs to use for the comparison of control components? Takes either a vector of pathway IDs, 
        #' or \code{all} to process all pathways, or \code{none} to not use this comparison at all.
        #' @param pwid_to_plot Which pathways to plot for subsequent visual inspection? Takes either a vector of pathway IDs, 
        #' or \code{same} to process the same pathways as for control component comparison, 
        #' or \code{test5} to plot 5 random pathways (can change to any number),
        #' or \code{none} to not use this comparison at all.
        #' @param new_biopax First BioPAX object.
        #' @param orig_biopax Second BioPAX object.
        #' @param new_biopax_tag First BioPAX Object tag.
        #' @param orig_biopax_tag BioPAX object.
        #' @param verbose Logical. Show all pertaining progress?
        #' @param pw_df_path 
        #' @param source_name BioPAX source name.
        
        #' @author 
        #' Ivan Grishagin
        
        st<-Sys.time()
        
        if(pwid_to_compare=="none"){
            pwid_to_compare<-
                NULL
        } else if(pwid_to_compare=="all"){
            pwid_to_compare<-
                read_excel(path = pw_df_path
                           ,col_types = rep("text",11)
                           ,sheet = 1) %>%
                filter(!is.na(biopax.Pathway.Name)
                       ,Source==source_name) %>%
                .$biopax.Pathway.ID
            
        }
        
        if(pwid_to_plot=="none"){
            pwid_to_plot<-
                NULL
        } else if(pwid_to_plot=="same"){
            pwid_to_plot<-
                pwid_to_compare
        } else if(length(grep("test"
                              ,pwid_to_plot))>0){
            pwnum<-
                gsub("test"
                     ,""
                     ,pwid_to_plot) %>%
                as.integer
            if(is.na(pwnum)){
                stop("MAIN_biopax_comparison: wrong number of test pathways!")
            }
            set.seed(123)
            pwid_to_plot<-
                load.biopax.pathways(new_biopax)$biopax.Pathway.ID %>%
                sample(pwnum)
            
            message("Plotting "
                    ,pwnum
                    ," random pathways.")
        }
        
        #plot biopax graphs
        if(!is.null(pwid_to_plot)){
                check_biopax(biopax=new_biopax
                                  ,pw_to_check=pwid_to_plot
                                  ,what="id"
                                  ,output_identifier=new_biopax_tag
                                  ,verbose = verbose)
                
                check_biopax(biopax=orig_biopax
                                  ,pw_to_check=pwid_to_plot
                                  ,what="id"
                                  ,output_identifier=orig_biopax_tag
                                  ,verbose = verbose)
        }
       
        #extract controllers and controlleds
        #and compare them
        if(is.null(pwid_to_compare)){
            return()
        }
        message("Comparing control components of "
                ,length(pwid_to_compare)
                ," pathways..."
                )
        status_vect<-
            pwid_to_compare %>%
            sapply(new_biopax=
                       new_biopax
                   ,orig_biopax=
                       orig_biopax
                   ,function(pwid
                             ,new_biopax
                             ,orig_biopax){
                       new_bp_contr<-
                           pathway2RegulatoryGraph_Rancho(biopax=new_biopax
                                                          ,pwid=pwid
                                                          ,returnGraph=FALSE
                           )
                       
                       orig_bp_contr<-
                           pathway2RegulatoryGraph_Rancho(biopax=orig_biopax
                                                          ,pwid=pwid
                                                          ,returnGraph=FALSE
                           )
                #if both don't have components (likely, something went wrong)
                #return "nocomponents"
                if(is.null(new_bp_contr) &
                   is.null(orig_bp_contr)){
                    return("nocomponents")
                }
                curr_status<-
                    identical(new_bp_contr
                              ,orig_bp_contr)
                return(curr_status)
            })
        
        #make and record the control component comparison status dataframe
        if(!is.null(status_vect)){
            status_df<-
                data_frame(pwid=pwid_to_compare
                           ,identical=status_vect)
            write.table(status_df
                        ,file = paste(Sys.Date()
                                      ,source_name
                                      ,"biopax comparison report.txt"
                        )
                        ,sep = "\t"
                        ,quote = FALSE
                        ,row.names = FALSE
                        ,col.names = TRUE)
        }
        et<-Sys.time()
        
        msg<-
            paste0(source_name
                   ," biopax comparison took "
                   ,round(difftime(et
                                   ,st
                                   ,units = "mins")
                          ,1)
                   ," minutes.")
        message(msg)
        write(msg
              ,"time_log.txt"
              ,sep = "\n"
              ,append=TRUE)
    }
######################################## MAIN_biopax_comparison ########################################
