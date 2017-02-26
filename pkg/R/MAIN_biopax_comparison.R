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
             ,groupnum=NULL
             ,ngroups=NULL
             
             
    ){
print("boo")
        #' @title
        #' MAIN -- Compare BioPAX Objects
        #' @description 
        #' Compare pathways in two BioPAX objects. 
        #' @details 
        #' Allows to compare all control components between two biopax objects.
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
        #' @param pw_df_path Path to the source table with pathway, id, name, and source values.
        #' @param source_name BioPAX public source name.
        #' @param groupnum Split pathway ids into \code{ngroups} groups. 
        #' Only process group #\code{groupnum} out of a total of \code{ngroups}.
        #' @param ngroups See \code{groupnum}.
        
        #' @author 
        #' Ivan Grishagin
        
        #establish initial key parameters, if not supplied by user
        if(pw_df_path=="default"){
            pw_df_path<-
                system.file("extdata"
                            ,"pathways_matched_to_sources_current_version.xlsx"
                            ,package="RIGbiopax")
        }
        
        st<-Sys.time()
        
        if(pwid_to_compare=="none"){
            pwid_to_compare<-
                NULL
        } else if(pwid_to_compare=="all"){
            pwid_to_compare<-
                read_excel_astext(path = pw_df_path
                                  #,col_types = rep("text",11)
                                  ,sheet = 1) %>%
                filter(!is.na(biopax.Pathway.Name)
                       ,Source==source_name) %>%
                .$biopax.Pathway.ID
            
        }
        #take only a chunk of pathways based on quantile split provided
        if (!is.null(groupnum)
            & !is.null(ngroups)){
            #get split intervals based on desired number of intervals 
            #and vector length
            quant_splits<-
                quantile(1:length(pwid_to_compare)
                         ,probs=seq(0,1
                                    ,1/ngroups)
                         ,type=1)
            #define interval span
            if(groupnum==1){
                start<-1
            } else {
                start<-
                    quant_splits[groupnum]+1
            }
            end<-
                quant_splits[groupnum+1]
            #cut out desired interval span
            pwid_to_compare<-
                pwid_to_compare[start:end]
            message("Comparing pathways #"
                    ,start
                    ," to #"
                    ,end
                    ," from "
                    ,source_name
                    ," biopax.")
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
                    return("no_control_components")
                }
                new_in_orig<-
                    round(100*sum(new_bp_contr %in% orig_bp_contr)/length(new_bp_contr) 
                          ,digits = 0) 
                    
                orig_in_new<-
                    round(100*sum(orig_bp_contr %in% new_bp_contr)/length(orig_bp_contr) 
                          ,digits = 0) 
                   
                curr_status<-
                    paste(new_in_orig
                          ,length(new_bp_contr)
                          ,orig_in_new
                          ,length(orig_bp_contr)
                          ,orig_bp_contr[!orig_bp_contr %in% new_bp_contr]
                          ,sep=" | ")
                
                return(curr_status)
            })
        
        #make and record the control component comparison status dataframe
        if(!is.null(status_vect)){
            
            status_df<-
                status_vect %>% 
                strsplit(split="\\|") %>% 
                do.call(rbind.data.frame
                        ,.) %>% 
                cbind.data.frame(pwid_to_compare
                                 ,.)
            colnames(status_df)<-
                c("pwid"
                  ,"new-in-orig, %"
                  ,"N_new_control_components"
                  ,"orig-in-new, %"
                  ,"N_orig_control_components"
                  ,"orig_not_in_new"
                  )

            write.table(status_df
                        ,file = paste(Sys.Date()
                                      ,source_name
                                      ,groupnum
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
