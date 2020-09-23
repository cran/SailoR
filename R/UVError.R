UVError <-
function(U,V,Ensembles=FALSE){
  ##UVError function
  ##INPUT
  #U:Reference
  #V:Model to test 
  #Ensembles:
  
  ##OUTPUT
  #$meanU: mean of the U dataset
  #$meanV: mean of the V dataset
  #$TotVarU: Total variance of the U dataset
  #$TotVarV: Total variance of the V dataset
  #$Eu: EOF matrix for the U dataset
  #$Ev: EOF matrix for the V dataset
  #$Rvu: Rotation matrix to express Ev as a rotation from Eu
  #$Sigmau: Matrix with standard deviations of U (2x2 matrix)
  #$Sigmav: Matrix with standard deviations of V (2x2 matrix)
  #$sdUx: Standard deviation of zonal component (U)
  #$sdUy: Standard deviation of meridional component (U)
  #$sdVx: Standard deviation of zonal component (V)
  #$sdVy: Standard deviation of meridional component (V)
  #$thetau: Angle of EOF1 (U) with zonal axis
  #$thetav: Angle of EOF1 (V) with meridional axis
  #$thetavu: Rotation angle from U EOFs to V EOFs
  #$RMSE: Rootmean square error between U and V
  #$Rvc2: vector correlation squared in 2D.(Breaker, Gemmill, Crossby, 1994, J. App. Met.)
  #doi: https://doi.org/10.1175/1520-0450(1994)033<1354:TAOATF>2.0.CO;2
  #following R vector correlation squared from CCA (Crossby et. al 1993, J. Atmos. Ocean. Tech. vol10)
  #$EccentricityU: eccentricity of the ellipses from the reference dataset
  #$EccentricityV: eccentricity of the ellipses from the model dataset
  #$congruenceEOF: Congruence coefficient (absolute value) for EOF1
  
  #If Ensembles = TRUE, n differents models (from mod)
  #are taken as a single model and the reference model is repeated n times.
  if(Ensembles==TRUE){
    nref<-nrow(U)
    nmod<-nrow(V)
    times2repeat<-nmod %/% nref
    
    if((nmod %% nref)!=0){
      message("--------------------------------------------")
      message("The number of elements in the ensemble is not
            an integer multiple of the number of observations,
            leaving")
      return(NULL)
    }
    
    U<-U[rep(1:nref,times2repeat),]
  }
  
  # Samples
  N=nrow(U)
  
  # Mean values
  meanU=colMeans(U)
  meanV=colMeans(V)
  
  # Reformat means to a full matrix
  muU=matrix(rep(meanU,N),nrow=N,byrow=TRUE)
  muV=matrix(rep(meanV,N),nrow=N,byrow=TRUE)
  
  # Anomalies
  Uanom=U-muU
  Vanom=V-muV
  
  # Principal components
  pcaU=prcomp(Uanom)
  pcaV=prcomp(Vanom)

  # Standard deviations of original data
  sdUx=sd(U[,1])
  sdUy=sd(U[,2])
  sdVx=sd(V[,1])
  sdVy=sd(V[,2])

  
  # Prepare decomposition matrices from PCA results
  Sigmau=matrix(rep(0,4),nrow=2)
  diag(Sigmau)=pcaU$sdev
  Sigmav=matrix(rep(0,4),nrow=2)
  diag(Sigmav)=pcaV$sdev
  TotVarU=sum(diag(Sigmau)*diag(Sigmau))
  TotVarV=sum(diag(Sigmav)*diag(Sigmav))
  # Eccentricities of ellipses
  EccentricityU=sqrt(1.-Sigmau[2,2]**2/Sigmau[1,1]**2)
  EccentricityV=sqrt(1.-Sigmav[2,2]**2/Sigmav[1,1]**2)

  
  # Principal components
  Pu=pcaU$x
  Pv=pcaV$x
  
  # E rotations
  Eu=pcaU$rotation
  Ev=pcaV$rotation
  # Angles with zonal axis
  thetau=atan2(Eu[2,1],Eu[1,1])
  # This is not needed There is a potential future problem with this
  # definition. In the equations we assume a relation between eof1 and eof2
  # so that eof2 is orthogonal to eof1 and we advance from eof1 to eof2
  # in the positive direction of the angle theta
  # BUT LAPACK routines do not ensure that, so better not to define
  # Ru this way
  # The right way would be
  # Eu=Ru%*%Identity, thus Ru=Eu, but ths is not needed, so, we remove it
  # Ru=matrix(c(cos(thetau),sin(thetau),-sin(thetau),cos(thetau)),nrow=2)
  thetav=atan2(Ev[2,1],Ev[1,1])
  # See above, unneeded but risky, do NOT use it or, if needed,
  # Rv=Ev
  # Rv=matrix(c(cos(thetav),sin(thetav),-sin(thetav),cos(thetav)),nrow=2)
  # Congruence coefficient
  e1u=Eu[,1]
  e1v=Ev[,1]
  congruenceEOF=abs(sum(e1u*e1v))

  # Rotation matrix.
  # This is a safe way of defining it (see comment in Ru above) and ensures
  # that equation (17) in the paper is true for any relative orientation of
  # eofs 1 and 2 
  Rvu=Ev%*%t(Eu)
  # This is also true if we focus on the relative orientation
  # between eof1(U) and eof1(V)
  thetavu=atan2(Rvu[2,1],Rvu[1,1])
  
  # Error matrix
  # Reworked using original data
  # Safe way to program the MSE, equation 21 in the paper
  DeltaUV<-function(A,B){
	N=nrow(A)
  	return(t(as.matrix(A-B))%*%as.matrix(A-B)/N)
  }
  
  # The frobenius norm of the error matrix in (21), so, eq (22)
  MSEError<-function(A,B){
	return(norm(DeltaUV(A,B),type="f"))
  }

  MSE=MSEError(U,V)
  
  # From MSE to RMSE
  RMSE=sqrt(MSE)
  
  #$Rvc2: R vector correlation squared from CCA 
  # (Crossby et. al 1993, J. Atmos. Ocean. Tech. vol10)
  # Canonical correlations following Wilks, Chapter 12, 12:30
  # Auxiliary function, matrix to 1/2
  Matrix2Half<-function(MM){
    KK=eigen(MM)
    retVal=KK$vectors%*%matrix(c(sqrt(KK$values[1]),0,0,sqrt(KK$values[2])),nrow=2,byrow=TRUE)%*%solve(KK$vectors)
    return(retVal)
  }
  
  # Inverse of matrix to 1/2 (matrix to -1/2)
  Matrix2InvHalf<-function(MM){
    return(solve(Matrix2Half(MM)))
  }
  
  # Compute CCA coefficients from data (anomalies computed inside)
  # Checked internal computations with example in wilks
  # Sxx=matrix(c(59.516,75.433,75.433,185.467),nrow=2,byrow=TRUE)
  # Syy=matrix(c(61.847,56.119,56.119,77.581),nrow=2,byrow=TRUE)
  # Syx=matrix(c(58.070,81.633,51.697,110.800),nrow=2,byrow=TRUE)
  # Wilks: 0.960 and 0.770
  CCAcoefficients<-function(UVr,UVm,mustCenter=TRUE){
    if (mustCenter){
      meanm=matrix(rep(colMeans(UVm),nrow(UVm)),ncol=2,byrow=TRUE)
      meanr=matrix(rep(colMeans(UVr),nrow(UVr)),ncol=2,byrow=TRUE)
      anomm=UVm-meanm
      anomr=UVr-meanr
    }else{
      anomm=UVm
      anomr=UVr
    }
    N=nrow(anomm)
    Sxx=t(anomm)%*%anomm/N
    Syy=t(anomr)%*%anomr/N
    Sxy=t(anomm)%*%anomr/N
    Syx=t(Sxy)
    lhs=Matrix2InvHalf(Sxx)%*%Sxy%*%Matrix2InvHalf(Syy)
    theSVD=svd(lhs)
    return(theSVD$d)
  }
  
  # This is the one called from the exterior
  # If anomalies are provided, mustCenter must be set to FALSE
  R2_Based_on_CCA<-function(UVr,UVm,mustCenter=TRUE){
    ccc=CCAcoefficients(UVr,UVm,mustCenter)
    return(sum(ccc**2))
  }
  
  R2vec<-R2_Based_on_CCA(as.matrix(Uanom),as.matrix(Vanom),FALSE)
  names(R2vec)<-NULL
  
  
  return(list(meanU=meanU,meanV=meanV,TotVarU=TotVarU,
              TotVarV=TotVarV,Eu=Eu,Ev=Ev,Rvu=Rvu,
              Sigmau=Sigmau,Sigmav=Sigmav,
              sdUx=sdUx,sdUy=sdUy,sdVx=sdVx,sdVy=sdVy,
              thetau=thetau,thetav=thetav,thetavu=thetavu,
              RMSE=RMSE,R2vec=R2vec,EccentricityU=EccentricityU,
	      EccentricityV=EccentricityV,congruenceEOF=congruenceEOF))
  }
