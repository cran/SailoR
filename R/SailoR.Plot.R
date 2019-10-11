SailoR.Plot <-
function(ref,mod,
                      ColourList=NULL,
                      sfactor=1,
                      docenter=FALSE,
                      xlim=NULL,
                      ylim=NULL,
                      xlab=NULL,
                      ylab=NULL,
                      plotmain=NULL,
                      plotlegend=TRUE,
                      plotRMSElegend=TRUE,
                      Ensembles=FALSE,
                      RMSE_legend_units="",
                      RMSE_legend_Rounding=2,
                      referenceName=NULL ){
  ##SailoR.plot function
  ##INPUT
  #ref: is data.frame with reference observations informationobject, 3 columns: model name, U component and V component
  #mod: is data.frame with model data, 3 columns: model name and U component, V component
  #ColourList: list with the colors you want to paint the models. Grey will always be used for reference data.
  #sfactor: 
  #docenter:
  #xlim:
  #ylim:
  #xlab:
  #ylab:
  #plotmain:
  #plotlegend: TRUE/ FALSE for automatical legend
  #plotRMSElegend: TRUE/FALSE plot another legend with RMSE results
  #Ensembles: TRUE/FALSE if Ensembles = TRUE, n differents models (from mod)
  #            are taken as a single model and the reference model is repeated n times.
  # RMSE_legend_units : Units to write in the RMSE legend
  # RMSE_legend_Rounding : Number of decimals to keep in RMSE legend
  # referenceName : External name provided for the reference dataset
  

  
  ##OUTPUT
  #$p: SailoR plot
  
  #PlotOneModel: plot the ellipses for one model
  PlotOneModel<-function(UV,thisColour,sfactor=1,docenter=FALSE,RMSE){
    # Plot meanV
    points(UV$meanV[1],UV$meanV[2],col=thisColour,pch=16,cex=1.1)
    
    #Plot RMSE
    add.alpha <- function(col, alpha=1){
      if(missing(col))
        stop("Please provide a vector of colours.")
      apply(sapply(col, col2rgb)/255, 2, 
            function(x) 
              rgb(x[1], x[2], x[3], alpha=alpha))  
    }
  
    
    # Plot the ellipse centered at this point
    # Prepare an ellipse in the space of the PCs
    nPoints=100
    # Phases
    dtheta=2*pi*1.01/nPoints
    theta=dtheta*(1:nPoints)
    
    # Ellipse in the space of the PCs (V) (axes)
    Pe=matrix(rep(c(0,0),nPoints),nrow=nPoints)
    Pe[,1]=UV$Sigmav[1,1]*cos(theta)*sfactor
    Pe[,2]=UV$Sigmav[2,2]*sin(theta)*sfactor
    # Transform back to physical space and center it into meanV
    Ptransf=Pe%*%t(UV$Ev)
    if(docenter==FALSE){
      Ptransf= Ptransf + matrix(rep(UV$meanV),nrow=nPoints,ncol=2,byrow=TRUE)
    }else{
      Ptransf= Ptransf + matrix(rep(UV$meanU),nrow=nPoints,ncol=2,byrow=TRUE)
    }
    lines(Ptransf,col=thisColour,lw=2)
    
    # Plot its major axis (as a fraction of variance of U!!, it is NOT a bug)
    Pe=matrix(rep(0,4),nrow=2)
    Pe[1,]=c(UV$Sigmav[1,1]*sfactor,0)
    Pe[2,]=c(-UV$Sigmav[1,1]*sfactor,0)
    Ptransf=Pe%*%t(UV$Ev)
    if(docenter==FALSE){
       Ptransf=Ptransf+matrix(rep(UV$meanV),nrow=2,ncol=2,byrow=TRUE)
     }else{
       Ptransf= Ptransf + matrix(rep(UV$meanU),nrow=2,ncol=2,byrow=TRUE)
    }
    lines(Ptransf,col=thisColour,lw=2)
    
    # Ellipse with the structure corresponding to U as reference (angular and axes diff)
    Pe=matrix(rep(c(0,0),nPoints),nrow=nPoints)
    Pe[,1]=UV$Sigmau[1,1]*cos(theta)*sfactor
    Pe[,2]=UV$Sigmau[2,2]*sin(theta)*sfactor
    Ptransf=Pe%*%t(UV$Eu)
    if(docenter==FALSE){
      Ptransf=Ptransf+matrix(rep(UV$meanV),nrow=nPoints,ncol=2,byrow=TRUE)
    }else{
      Ptransf= Ptransf + matrix(rep(UV$meanU),nrow=nPoints,ncol=2,byrow=TRUE)
    }
    lines(Ptransf,col="gray",lw=2,type="l",lt=2)
    
    # Major axis of this ellipse
    Pe=matrix(rep(0,4),nrow=2)
    Pe[1,]=c(UV$Sigmau[1,1]*sfactor,0)
    Pe[2,]=c(-UV$Sigmau[1,1]*sfactor,0)
    Ptransf=Pe%*%t(UV$Eu)
    if(docenter==FALSE){
      Ptransf=Ptransf+matrix(rep(UV$meanV),nrow=2,ncol=2,byrow=TRUE)
    }else{
      Ptransf= Ptransf + matrix(rep(UV$meanU),nrow=2,ncol=2,byrow=TRUE)
    }
    lines(Ptransf,col="gray",lw=2,lt=2)
  }
  
  #If Ensembles = TRUE, n differents models (from mod)
  #are taken as a single model and the reference model is repeated n times. 
  # if(Ensembles==TRUE){
  #   ref<-ref[rep(1:nrow(ref),length(unique(mod$mod))),]
  #   mod$mod<-"Ens"
  # }
  
  UV<-SailoR.Indices(ref,mod,Ensembles)
  # Plot means and axes
  plot(UV[[1]]$meanU[1],UV[[1]]$meanU[2],type="p",
        xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,
        pch=16,col="darkgray",
        cex=1.1,main=plotmain,asp=1)
  grid()
    
  if(is.null(ColourList)==TRUE){
    ColourList<-c("darkgrey","red","blue","darkolivegreen3","orange","seagreen3","gold","purple","pink")
    #ColourList<-c("blue","red",'springgreen4',"orange","gray50","violetred","orchid2","chocolate")
  }else{
     ColourList<-c("darkgrey",ColourList)
  }
    
  #RMSE for all models
  RMSE<-NULL
  for(i in 1:length(UV)){
    RMSE<-c(RMSE,UV[[i]]$RMSE)
  }
    
  for(i in 1:length(UV)){
    PlotOneModel(UV[[i]],ColourList[i+1],sfactor,docenter,RMSE)
  }
    
  # Plot meanU (again, on top of everything, just in case)
  points(UV[[1]]$meanU[1],UV[[1]]$meanU[2],col="darkgray",pch=15,cex=1.1)
    
  # Legend
  if(plotlegend==TRUE){
    
    
    ############## AQUÍ
    if(is.null(referenceName)){
      lrefname=as.character(ref$mod[1])
    }else{
      lrefname=referenceName
    }
    
    legend(x="topright",
        c(lrefname,names(UV)),
        col=c(ColourList[c(1:(length(UV)+1))]),lwd=2)
  }
  
  #RMSE legend
  if(plotRMSElegend==TRUE){
    legend(x="bottomright",
           paste(names(UV),": ",round(RMSE,RMSE_legend_Rounding),RMSE_legend_units,sep=""),
           col=ColourList[c(2:(length(UV)+1))],
           pch=19,title=expression(bold("RMSE")))  
  }
  
  #plot
  p <- recordPlot()
  return (p)
}
