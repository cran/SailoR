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
                      plotscalelegend=TRUE,
                      Ensembles=FALSE,
                      RMSE_legend_units="",
                      RMSE_legend_Rounding=2,
                      referenceName=NULL,
                      linestype=NULL,
                      bias_pch=NULL){
  ##SailoR.plot function
  ##INPUT
  #ref: is data.frame with reference observations informationobject, 3 columns: model name, U component and V component
  #mod: is data.frame with model data, 3 columns: model name and U component, V component
  #ColourList: list with the colors you want to paint the models. Dark gray will always be used for reference data.
  #The default option provides 10 colors of virids palette which is robust to colorblindness.
  #sfactor: 
  #docenter:
  #xlim:
  #ylim:
  #xlab:
  #ylab:
  #plotmain:
  #plotlegend: TRUE/ FALSE for automatical legend
  #plotRMSElegend: TRUE/FALSE plot another legend with RMSE results
  #plotscalelegend: TRUE/FALSE plot another legend wiht scaled factor
  #Ensembles: TRUE/FALSE if Ensembles = TRUE, n differents models (from mod)
  #            are taken as a single model and the reference model is repeated n times.
  # RMSE_legend_units : Units to write in the RMSE legend
  # RMSE_legend_Rounding : Number of decimals to keep in RMSE legend
  # referenceName : External name provided for the reference dataset
  #bias_pch: a vector with the pch types for the bias simbol.
  # The default opcion provides a circle when docenter=TRUE and point when docenter=FALSE

  
  ##OUTPUT
  #$p: SailoR plot
  
  #PlotOneModel: plot the ellipses for one model
  PlotOneModel<-function(UV,thisColour,thislinestype,sfactor=1,docenter=FALSE,RMSE,POINTS,thisbias_pch){
    # Plot meanV
    if(docenter==TRUE){
        cex=3.1
        pch=thisbias_pch
    }else{
        cex=2
        if(thisbias_pch==1){
            pch=16
        }else{pch=thisbias_pch}
    }
    if(POINTS==TRUE){
        points(UV$meanV[1],UV$meanV[2],col=thisColour,pch=pch,lwd=3,cex=cex)
    }else{

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
        lines(Ptransf,col=thisColour,lwd=2,lty=thislinestype)
        
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
        lines(Ptransf,col=thisColour,lwd=2,lty=thislinestype)
        
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
        lines(Ptransf,col="gray40",lwd=2,type="l",lty=2)
        
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
        lines(Ptransf,col="gray40",lwd=2,lty=2)
    }
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
        pch=1,lwd=3,col=NA,
        cex=2,main=plotmain,asp=1)
  grid()
    
  if(is.null(ColourList)==TRUE){
    ColourList<-c("gray40","red","blue","darkolivegreen3",
                  "orange","seagreen3","gold","purple","pink")

    # ColourList<-c("gray40","#FDE725FF","#482878FF",
    #                         "#6DCD59FF","#31688EFF",
    #                         "#B4DE2CFF","#3E4A89FF",
    #                         "#35B779FF","#26828EFF",
    #                         "#1F9E89FF")
  }else{
     ColourList<-c("gray40",ColourList)
  }
  
  
 
  n<-which(linestype==2 | linestype=="dashed")
  if(length(n)!=0){
    linestype[n]<-"solid"
  }
  
  if(is.null(linestype)==FALSE){
      if(length(linestype)!=length(UV)){
          linestype<-rep(linestype,length.out=length(UV))
      }
  }else{
      linestype<-rep("solid",length.out=length(UV))
  }
  
  if(is.null(bias_pch)==FALSE){
    if(length(bias_pch)!=length(UV)){
      bias_pch<-rep(bias_pch,length.out=length(UV))
    }
  }else{
    bias_pch<-rep(1,length.out=length(UV))
  }
  
  #RMSE for all models
  RMSE<-NULL
  for(i in 1:length(UV)){
    RMSE<-c(RMSE,UV[[i]]$RMSE)
  }
    
  for(i in 1:length(UV)){
    PlotOneModel(UV[[i]],ColourList[i+1],linestype[i],sfactor,docenter,RMSE,POINTS=FALSE,bias_pch[i])
  }
  
  for(i in 1:length(UV)){
    PlotOneModel(UV[[i]],ColourList[i+1],linestype[i],sfactor,docenter,RMSE,POINTS=TRUE,bias_pch[i])
  }  
  
  # Plot meanU (again, on top of everything, just in case)
  points(UV[[1]]$meanU[1],UV[[1]]$meanU[2],col="gray40",pch=15,cex=2)
    
  # Legend
  if(plotlegend==TRUE){
    
    if(is.null(referenceName)){
      lrefname=as.character(ref$mod[1])
    }else{
      lrefname=referenceName
    }
    
    legend(x="topright",
        c(lrefname,names(UV)),
        col=c(ColourList[c(1:(length(UV)+1))]),lwd=2,
        lty=c("dashed",linestype),cex=1,ncol=1)
  }
  
  #RMSE legend
  if(plotRMSElegend==TRUE){
    legend(x="bottomright",
           paste(names(UV),": ",round(RMSE,RMSE_legend_Rounding),
                 RMSE_legend_units,sep=""),
           col=ColourList[c(2:(length(UV)+1))], lty=linestype,
           pch=bias_pch,lwd=3,title=expression(bold("RMSE")),
           cex=1,ncol=1) 
  }
  

  #scale legend
  if(plotscalelegend==TRUE){
    par(font=2)
    legend(x="topleft",paste("Scaled with ",sfactor,sep=""), bty = "n",ncol=1,cex=1) 
  }
  
  
  #plot
  p <- recordPlot()
  return (p)
}
