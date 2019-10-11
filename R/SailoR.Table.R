SailoR.Table <-
function(UV,round_digits=2){
  
    ##Index.Table function
    ##INPUT
    #UV: SailoR.Indices output
  
    ##OUTPUT
    #Table: is a data.frame with the SailoR.Indices information in table format
  
    #Function
    euclideanNorm<-function(a){
      return(sqrt(a[1]*a[1]+a[2]*a[2]))
    }
    
    keyWords<-c("Ref",names(UV))
    sdUx<-UV[[names(UV)[1]]]$sdUx
    sdUy<-UV[[names(UV)[1]]]$sdUy
    sdVx<-NA
    sdVy<-NA
    thetau<-UV[[names(UV)[1]]]$thetau
    thetav<-NA
    thetavu<-NA
    R2vec<-2
    biasMag<-0
    RMSE<-0
    
    for(i in 1:length(names(UV))){
        sdUx<-c(sdUx,NA)
        sdUy<-c(sdUy,NA)
        sdVx<-c(sdVx,UV[[names(UV)[i]]]$sdVx)
        sdVy<-c(sdVy,UV[[names(UV)[i]]]$sdVy)
        thetau<-c(thetau,NA)
        thetav<-c(thetav,UV[[names(UV)[i]]]$thetav)
        thetavu<-c(thetavu,UV[[names(UV)[i]]]$thetavu)
        R2vec<-c(R2vec,UV[[names(UV)[i]]]$R2vec)
        biasMag<-c(biasMag,euclideanNorm(UV[[names(UV)[i]]]$meanU-UV[[names(UV)[i]]]$meanV))
        RMSE<-c(RMSE,UV[[names(UV)[i]]]$RMSE)
    }
    
    TABLE<-data.frame(keyWords,round(sdUx,round_digits),round(sdUy,round_digits),
                      round(sdVx,round_digits),round(sdVy,round_digits),
                      round(thetau,round_digits),round(thetav,round_digits),
                      round(thetavu,round_digits),
                      round(R2vec,round_digits),
                      round(biasMag,round_digits),round(RMSE,round_digits))
    colnames(TABLE)<-c("modelName","sdUx","sdUy","sdVx","sdVy",
                       "thetau","thetav","thetavu",
                       "R2vec","biasMag","RMSE")
    return(TABLE)
}
