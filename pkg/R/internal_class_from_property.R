internal_class_from_property<-
    function(property){
        #returns a class based on the supplied property
        #some cases are very well-defined, so this is justified
        switch(property
               ,component = "Complex"
               ,modificationType = "ModificationFeature"
               ,featureLocation = "FragmentFeature"
               ,evidenceCode = "Evidence"
               ,sequenceIntervalBegin = "SequenceInterval"
               ,sequenceIntervalEnd = "SequenceInterval"
               ,"PLACEHOLDER")
    }