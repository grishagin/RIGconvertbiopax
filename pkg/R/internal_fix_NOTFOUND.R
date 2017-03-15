internal_fix_NOTFOUND<-
    function(biopax_dt){
        if(! "NOTFOUND" %in% biopax_dt$class){
            return(biopax_dt)
        } 
        biopax_dt<-
            as.data.table(biopax_dt)
        
        notfound_ids<-
            biopax_dt[class %in% c("NOTFOUND"
                                   ,"NOTFOUNDReference")]$id
        
        message("internal_fix_NOTFOUND: Found 'NOTFOUND' class items with ids "
                ,paste(notfound_ids
                       ,collapse=", ")
                ,"\nTrying to fix known ones...")
        
        ########## fix keNOTFOUND1 -- missing protein name in KEGG
        biopax_dt$class[biopax_dt$id=="keNOTFOUND1"]<-
            "Protein"
        biopax_dt$class[biopax_dt$id=="keNOTFOUNDReference1"]<-
            "ProteinReference"
        
        new_entries<-
            data.table(class="Protein"
                       ,id="keNOTFOUND1"
                       ,property="displayName"
                       ,property_attr="rdf:resource"
                       ,property_attr_value="http://www.w3.org/2001/XMLSchema#string"
                       ,property_value="MHC class I polypeptide-related sequence A")
        
        biopax_dt<-
            rbind.data.frame(biopax_dt
                             ,new_entries)
        
        ########## check if there are more
        if("NOTFOUND" %in% biopax_dt$class){
            msg<-
                "internal_fix_NOTFOUND: More 'NOTFOUND' class items found, you should check the final output!"
            message(msg)
            warning(msg)
            
        } 
        return(biopax_dt)
    }