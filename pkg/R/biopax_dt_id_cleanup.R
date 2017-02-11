######################################## biopax_dt_id_cleanup ########################################
biopax_dt_id_cleanup<-
    function(biopax_dt){
        
        #' @title
        #' Clean BioPAX Data Table
        #' @description 
        #' Takes in a BioPAX-style datatable, 
        #' and cleans it by removing duplicate entities, that have the same values, but different ids. 
        #' @details 
        #' Outputs a BioPAX-style datatable. 
        #' @param biopax_dt BioPAX-style datatable.
        
        #' @author 
        #' Ivan Grishagin
        
        
        #function to remove duplicate entities
        #that have the same values, but different ids
        
        #first, select entries with property value not equal to blank
        #and add a column of combined values
        #for ease of processing
        allowed_props<-
            c("name"
              ,"displayName"
              ,"standardName"
              ,"id")
        
        dT<-
            biopax_dt %>%
            filter(nchar(property_value)>0
                   ,property %in% allowed_props) %>%
            mutate(temp=paste0(class
                               ,property
                               ,property_attr
                               ,property_attr_value
                               ,property_value))
        #entries that have duplicates
        dupl_entr<-
            dT$temp %>%
            table(.) %>%
            .[.>1] %>%
            names
        
        #repace values
        for(entr in dupl_entr){
            to_repl<-
                dT %>%
                filter(temp %in% entr) %>%
                .$id
            biopax_dt$id[biopax_dt$id %in% to_repl]<-
                to_repl[1]
            biopax_dt$property_attr_value[biopax_dt$property_attr_value %in% to_repl]<-
                to_repl[1]
        } 
        biopax_dt<-
            biopax_dt %>%
            unique
        
        return(biopax_dt)
        
    }
######################################## biopax_dt_id_cleanup ########################################
