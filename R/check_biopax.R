######################################## check_biopax ########################################
check_biopax<-
    function(biopax
             ,pw_to_check=NULL
             ,what=c("sample","order","id")
             ,output_identifier=NULL
             ,verbose=TRUE){
        
        #' @title
        #' Check BioPAX by Plotting Certain Pathways
        #' @description 
        #' Plot regulatory graph for each supplied pathway.
        #' @param biopax BioPAX object.
        #' @param pw_to_check Which pathways to check? Can take in (1) pathway ids or
        #' (2) integer as number of pathways to randomly select and check 
        #' (3) vector of integers as pathway indices.
        #' @param what What's the format of supplied ids? Either \code{id}, or \code{sample}, or \code{order}.
        #' @param output_identifier What to append to output file?
        #' @param verbose Logical. Show all pertaining progress?
        
        #' @author 
        #' Ivan Grishagin
        
        require(graph)
        require(Rgraphviz)
        if(is.null(pw_to_check)){
            message("check_biopax: pw_to_check=NULL. Won't plot anything.")
            return()
        }
        
        if("character" %in% class(biopax)){
            biopax<-
                readBiopax(biopax)
        } else if (!("biopax" %in% class(biopax))) {
            stop("check_biopax: parameter biopax is neither a filename nor a biopax object.")
        }
        all_pwid<-
            listPathways(biopax)$id
        
        if(what=="sample"){
            pwids<-
                sample(all_pwid
                       ,pw_to_check)
        } else if(what=="order"){
            if(any(pw_to_check>length(all_pwid)) | pw_to_check<1){
                stop("check_biopax: order vector pw_to_check = c("
                     ,pw_to_check
                     ,"), while there are "
                     ,length(all_pwid)
                     ," pathways in biopax.")
            }
            pwids<-
                listPathways(biopax)$id[pw_to_check]
        } else if(what=="id"){
            if(any(!(pw_to_check %in% all_pwid))){
                stop("check_biopax: id vector pw_to_check = c("
                     ,pw_to_check
                     ,") has ids that are not present in biopax.")
            }
            pwids<-
                pw_to_check
        } 
        
        for (pwid in pwids){
            message("Making and plotting pathway "
                    ,pwid
                    ,output_identifier
                    ," graph...")
            graph<-
                pathway2RegulatoryGraph_Rancho(biopax = biopax
                                               ,pwid = pwid
                                               ,verbose = verbose
                                               ,returnGraph=TRUE)
            
            if(is.null(graph)){
                message("check_biopax: pwid "
                        ,pwid
                        ," does not seem to have any components to plot.")
                next()
            }
            
            png(filename = paste0(Sys.Date()
                                  ,"_biopax_graph_pw_"
                                  ,pwid
                                  ,"_",output_identifier
                                  ,".png")
                ,width = 8400
                ,height = 8400
                ,res=2400
                ,pointsize = 15)

            tryCatch(plotRegulatoryGraph(graph)
                     ,error = function(err) {
                         message("plotRegulatoryGraph: error in pathway "
                                 ,pwid
                                 ,".\n"
                                 ,err)
                     })
            dev.off()
            
        }
    }
######################################## check_biopax ########################################
