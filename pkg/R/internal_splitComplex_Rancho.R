internal_splitComplex_Rancho<-
    function (biopax
              ,complexid
              ,recursive = TRUE
              ,returnIDonly = FALSE
              ,biopaxlevel = 3) {
        if ("biopax" %in% class(biopax)) {
            df<-
                biopax$dt
            biopaxlevel<-
                biopax$biopaxlevel
        }
        else if ("data.table" %in% class(biopax)) {
            df<-
                biopax
        }
        else {
            stop("splitComplex: parameter biopax is neither biopax object nor compatible biopax data.table")
        }
        compname<-
            c("COMPONENTS"
              ,"PHYSICAL-ENTITY")
        if (biopaxlevel == 3) {
            compname<-c("component"
                        #added by Rancho
                        ,"memberPhysicalEntity")
        }
        ref<-
            getReferencedIDs(df
                             ,complexid
                             ,recursive = recursive
                             ,onlyFollowProperties = compname)
        if (is.null(ref)) 
            return(NULL)
        referenced<-
            selectInstances(df
                            ,id = ref
                            ,returnCopy = FALSE
                            ,biopaxlevel = biopaxlevel) %>% 
            filter(tolower(class) %chin% 
                       c("dna"
                         ,"rna"
                         ,"protein"
                         ,"smallmolecule")
                   ,tolower(property) %chin%
                       c("name"
                         ,"displayname"
                         ,"standardname"
                         ,"term")) %>% 
            dplyr::select(id
                          ,name=property_value)
        # referenced<-
        #     unique(as.character(referenced[tolower(class) %chin% 
        #                                        c("dna"
        #                                          ,"rna"
        #                                          ,"protein"
        #                                          ,"smallmolecule")]$id))
        if (length(referenced) == 0) {
            return(NULL)
        }
        if (returnIDonly) {
            return(striphash(referenced$id))
        }
        # return(listInstances(df
        #                      ,id = referenced
        #                      ,biopaxlevel = biopaxlevel))
        return(referenced)
    }