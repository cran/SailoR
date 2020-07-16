SailoR.Indices <-
function(ref,mod,Ensembles=FALSE){

    ##INPUT
    #ref: is data.frame with reference observations informationobject, 3 columns: model name, U component and V component
    #mod: is data.frame with model data, 3 columns: model name and U component, V component
    #Ensembles: TRUE if model includes ensemble-runs
  
    ##OUTPUT
    #$result: a data.frame with columns "meanU","meanV","RMSE",
    #"sdUx","sdUy","sdVx","sdVy","Eu","Ev",
    #"Sigmau","Sigmav","TotVarU","TotVarV","thetau","thetav",
    #"thetavu","R2vec","Rvu","EccentricityU",
    #"EccentricityV","congruenceEOF"
    #With names as described in UVError() output
  
    #Transformation data and check col dimension if Ensembles is FALSE
    nref<-nrow(ref)
    nmod<-nrow(mod)
    if((nmod %% nref)!=0){
      message("--------------------------------------------")
      message("The number of elements in the models is not
            an integer multiple of the number of observations,
            leaving")
      return(NULL)
    }
    
    #Check U and V data class
    if(class(ref[,2])!="numeric" | class(ref[,3])!="numeric"){
        message("U, V columns in the reference have to be 'numeric' class")
        return(NULL)
    }
    
    if(class(mod[,2])!="numeric" | class(mod[,3])!="numeric"){
      message("U, V columns in the reference have to be 'numeric' class")
      return(NULL)
    }

    
    if(is.null(mod)==FALSE){
      
      #Calculation of UVError  
      UV<-list()
      n_jj<-unique(mod$mod)
      for(jj in 1:length(n_jj)){
        n_mod<-which(mod$mod==n_jj[jj])
        uverror<-UVError(ref[,c(2,3)],mod[n_mod,c(2,3)],Ensembles)
        if(is.null(uverror)==TRUE){
            message("Error in UVError")
            return(NULL)
        }else{
            UV[[jj]]<-uverror
            names(UV)[jj]<-as.character(n_jj[jj])
        }
      }
      
      #List with final object
      l<-list()
      for(ii in 1:length(UV)){
        l[[ii]]<-UV[[ii]][c("meanU","meanV","RMSE","sdUx","sdUy","sdVx","sdVy","Eu","Ev",
                            "Sigmau","Sigmav","TotVarU","TotVarV","thetau","thetav",
                            "thetavu","R2vec","Rvu","EccentricityU",
                            "EccentricityV","congruenceEOF")]
      }
      names(l)<-names(UV)
      return (l)
    }else{
      message("The dimensions don't match")
    }
}
