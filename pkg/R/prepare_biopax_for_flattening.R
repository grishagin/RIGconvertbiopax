prepare_biopax_for_flattening<-
    function(biopax){
        
        #' @title
        #' Prepare BioPAX for Flattening
        #' @description 
        #' Prepares a BioPAX object for correct flattening:  \cr
        #' 1. Removes hashes from property_attr_value column;  \cr
        #' 2. For instances of class \code{Complex}, replaces \code{memberPhysicalEntity} with \code{component};  \cr
        #' 3. Removes all \code{memberPhysicalEntity}-related instances by replacing 
        #' references to instances with \code{memberPhysicalEntity} property with references to the instances they reference.
        #' @param biopax BioPAX object.
        #' 
        #' @author 
        #' Ivan Grishagin
        
        #strip hashes from biopax property attr value column
        biopax$dt$property_attr_value<-
            biopax$dt$property_attr_value %>% 
            striphash
        
        #replace memberPhysicalEntity with complex component
        #for each such case, replace memberPhysicalEntity with component -- 
        #because that's what it essentially is
        mem_phys_ent_complex_logi<-
            biopax$dt$property=="memberPhysicalEntity" &
            biopax$dt$class=="Complex"
        biopax$dt$property[mem_phys_ent_complex_logi]<-
            "component"
        
        biopax<-
            biopax %>% 
            internal_remove_biopax_property(property_to_remove="memberPhysicalEntity")
        
        return(biopax)
    }