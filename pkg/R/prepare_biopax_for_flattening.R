prepare_biopax_for_flattening<-
    function(biopax){
        
        #' @title
        #' Prepare BioPAX for Flattening
        #' @description 
        #' Prepares a BioPAX object for correct flattening:  \cr
        #' 1. Applies \code{clean_biopax} function;  \cr
        #' 2. For instances of class \code{Complex}, replaces \code{memberPhysicalEntity} with \code{component};  \cr
        #' 3. Removes all \code{memberPhysicalEntity}-related instances by replacing 
        #' references to instances with \code{memberPhysicalEntity} property with references to the instances they reference.
        #' @param biopax BioPAX object.
        #' 
        #' @author 
        #' Ivan Grishagin
        
        #replace memberPhysicalEntity with complex component
        #for each such case, replace memberPhysicalEntity with component -- 
        #because that's what it essentially is
        biopax$dt[property=="memberPhysicalEntity" &
                      class=="Complex"]$property<-
            "component"
        
        biopax<-
            biopax %>% 
            #do some SERIOUS biopax cleaning
            clean_biopax %>% 
            #replace/reconnect memberPhysicalEntity
            internal_remove_biopax_property(property_to_remove="memberPhysicalEntity")
            
        
        return(biopax)
    }