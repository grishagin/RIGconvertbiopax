######################################## extract_biopax_interactions_vect ########################################
extract_biopax_interactions_vect<-
    function(pw_id_vect="all"
             ,biopax
             ,biopax_source
             ,tag=NULL){
        
        #' @title
        #' Extract BioPAX Interactions
        #' @description 
        #' Takes in a BioPAX object and returns a dataframe of all "interactions", one per row, for a specified vector of pathways.
        #' @details 
        #' Also outputs the resultant dataframe as *.txt and *.xlsx files. 
        #' @param pw_id_vect Vector of pathway ids to process. Can also take a value \code{all} to process all pathways.
        #' @param biopax BioPAX object.
        #' @param biopax_source Name of the public source of the BioPAX object.
        #' @param tag String to append to the output file.
        
        #' @author 
        #' Ivan Grishagin
        
        #function to extract interactions from biopax files
        #currently supported:
        #Biocarta
        #KEGG
        #NCI
        #Reactome
        #NetPath
        #Wiki
        require("dplyr")
        require("openxlsx")
        require("xlsx")
        require("rBiopaxParser")
        
        if(length(pw_id_vect)==1){
            if (pw_id_vect=="all"){
                pw_id_vect<-
                    listPathways(biopax = biopax
                                 ,biopaxlevel = biopax$biopaxlevel)$id
            }
        }
        
        #for each pathway, find all components, and 
        #extract interactions for them
        result_df_list<-
            lapply(pw_id_vect
                   ,extract_biopax_interactions
                   ,biopax=biopax
                   ,biopax_source=biopax_source)
        result_df<-
            do.call(rbind,result_df_list)
        
        try(write.table(result_df
                        ,paste(Sys.Date()
                               ,biopax_source
                               ,tag
                               ,"interactions.txt"
                               ,sep="_")
                        ,quote=FALSE
                        ,sep = "\t"
                        ,row.names = FALSE
                        ,col.names = TRUE))
        
        try(openxlsx::write.xlsx(result_df
                                 ,paste(Sys.Date()
                                        ,biopax_source
                                        ,tag
                                        ,"interactions.xlsx"
                                        ,sep="_")
                                 ,row.names = FALSE
                                 ,col.names = TRUE))
        
        return(result_df)
    }
######################################## extract_biopax_interactions_vect ########################################
