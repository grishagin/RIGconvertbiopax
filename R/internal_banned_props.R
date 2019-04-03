######################################## internal_banned_props ########################################
internal_banned_props<-
    function(){
        #banned properties (exclude from consideration)
        return(c("name"
                 ,"entityReference"
                 ,"xref"
                 ,"feature"
                 ,"modificationType"
                 ,"featureLocation"
                 ,"sequenceIntervalBegin"
                 ,"sequenceIntervalEnd"
                 ,"cellularLocation")
        )
    }
######################################## internal_banned_props ########################################
