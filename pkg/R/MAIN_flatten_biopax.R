######################################## MAIN_flatten_biopax ########################################
MAIN_flatten_biopax<-
    function(pwid_to_convert=c("test5"
                               ,"all")
             ,work_dir="default"
             ,pw_df_path="default"
             ,biopax=NULL
             ,source_name=NULL
             ,groupnum=NULL
             ,ngroups=NULL
    ){
        #' @title
        #' MAIN -- Flatten BioPAX 
        #' @description 
        #' Flatten BioPAX ontology file by extracting all "interactions". 
        #' @details 
        #' Returns a dataframe with one complete "interaction" (pathway component) per row. 
        #' Each "interaction" contains over 120 columns, flattening the BioPAX ontology as much as it is practically feasible. 
        #' Each complex component is represented as a special string structure, similar to JSON.
        #' @param pwid_to_convert Which pathway IDs to process? Takes either a vector of pathway IDs, 
        #' or \code{test5} to process 5 random pathways (can change to any number), 
        #' or \code{all} to process all pathways.
        #' @param pw_df_path Path to the source table with pathway, id, name, and source values.
        #' @param biopax BioPAX object.
        #' @param source_name BioPAX public source name.
        #' @param groupnum Split pathway ids into \code{ngroups} groups. 
        #' Only process group #\code{groupnum} out of a total of \code{ngroups}.
        #' @param ngroups See \code{groupnum}.
        
        #' @author 
        #' Ivan Grishagin
        
        #establish key parameters, if not supplied by user
        if(work_dir=="default"){
            work_dir<-getwd()
        }
        if(pw_df_path=="default"){
            pw_df_path<-
                system.file("extdata"
                            ,"pathways_matched_to_sources_current_version.xlsx"
                            ,package="RIGbiopax")
        }
        
        ########################## prepare
        if (!require("RIGessentials")){
            if(!require("devtools")){
                install.packages("devtools")
            }
            if(!require("devtools")){
                stop("Tried installing devtools and something is not working!")
            } 
            install_github("grishagin/RIGessentials"
                           ,subdir="pkg")
            if (!require("RIGessentials")){
                stop("Tried installing RIGessentials and something is not working!")
            }
        }
        
        pkg<-c("rJava"
               ,"plyr"
               ,"dplyr"
               ,"rBiopaxParser"
               ,"readxl"
               ,"openxlsx"
               ,"qdap"
               ,"stringr"
               ,"graph"
               ,"Rgraphviz")
        loadPackages(pkg
                     ,verbose = FALSE)
        
        if (!require("RIGbiopax")){
            if(!require("devtools")){
                install.packages("devtools")
            }
            if(!require("devtools")){
                stop("Tried installing devtools and something is not working!")
            } 
            install_github("grishagin/RIGbiopax"
                           ,subdir="pkg")
            if (!require("RIGbiopax")){
                stop("Tried installing RIGbiopax and something is not working!")
            }
        }
        
        prepareSession(work_dir
                       ,nolocale = FALSE)
        
        if(is.null(source_name)){
            #biopax sources
            biopax_source_names<-
                c("BioCarta"
                  ,"KEGG"
                  ,"NCI-Nature"
                  ,"NetPath"
                  ,"Wiki Pathways"
                  ,"Evgeny"
                  ,"Ayesha"
                  ,"Science Signaling"
                  ,"Reactome"
                )
            
            #source name
            source_name<-
                tkradio_from_vect(biopax_source_names
                                  ,"Select BioPAX Source.")
        }
        message("Processing "
                ,source_name
                ," biopax...")
        
        ########################## prepare 
        ########################## get pathways and biopax
        #select only pertaining pathways
        pw_df<-
            read_excel(path = pw_df_path
                       ,col_types = rep("text",11)
                       ,sheet = 1) %>%
            filter(!is.na(biopax.Pathway.Name)
                   ,Source==source_name) %>%
            dplyr::select(toxdb.Pathway.ID
                          ,toxdb.Pathway.Name
                          ,biopax.Pathway.Name
                          ,biopax.Pathway.ID
                          ,Source)
        
        #check if biopax object exists in environment, 
        #and if not -- load biopax
        if(is.null(biopax)){
            biopax<-
                load.biopax(source_name=source_name
                            ,source_dir=NULL)
        }
        
        ########################## get pathways and biopax
        ########################## clean biopax
        biopax<-
            clean_biopax_property_value(biopax)
        ########################## clean biopax
        ########################## pathways to convert
        if(length(grep("test"
                       ,pwid_to_convert))>0){
            pwnum<-
                gsub("test"
                     ,""
                     ,pwid_to_convert) %>%
                as.integer
            if(is.na(pwnum)){
                stop("MAIN_flatten_biopax: wrong number of test pathways!")
            }
            pwid_to_convert<-
                sample(pw_df$biopax.Pathway.ID
                       ,pwnum)
            message("Extracting "
                    ,pwnum
                    ," random pathways.")
        } else if(any(pwid_to_convert %in% "all")){
            pwid_to_convert<-
                pw_df$biopax.Pathway.ID
            message("Extracting all "
                    ,length(pwid_to_convert)
                    ," inxight pathways from "
                    ,source_name
                    ," biopax.")
        }
        #take only a chunk of pathways based on quantile split provided
        if (!is.null(groupnum) &
            !is.null(ngroups)){
            #get split intervals based on desired number of intervals 
            #and vector length
            quant_splits<-
                quantile(1:length(pwid_to_convert)
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
            pwid_to_convert<-
                pwid_to_convert[start:end]
            message("Extracting pathways #"
                    ,start
                    ," to #"
                    ,end
                    ," from "
                    ,source_name
                    ," biopax.")
        }
        ########################## pathways to convert
        ########################## extract interactions from biopax
        st<-Sys.time()
        interactions_df<-
            extractBiopaxInteractions_vect(pw_id_vect=pwid_to_convert
                                           ,biopax=biopax
                                           ,biopax_source=source_name
                                           ,tag=groupnum)
        
        et<-Sys.time()
        msg<-
            paste0("Conversion of "
                   ,length(pwid_to_convert)
                   ," "
                   ,source_name
                   ," pathways to flat file took "
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
       
        ########################## extract interactions from biopax
        try(invisible(tkmessageBox(message = "Conversion done!"
                                   ,icon = "info"
                                   ,type = "ok")))
        return(interactions_df)
    }
######################################## MAIN_flatten_biopax ########################################
