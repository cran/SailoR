\name{Ensembles}
\alias{Ensembles}
\docType{data}
\title{
Multimodel Ensembles
}
\description{
Contains the information of both components of surface wind over the Southern Hemisphere during the period 1979-2015 from the historical forcing experiment of CMIP5 repository. Four ensemble members from the IPSL model and three from the MIROC model are included. Additionally, data from the reanalysis ERA5 is also included. All of them have been bilinearly interpolated to a common 1.25º x 1º regular grid. 
}
\usage{data("Ensembles")}
\format{
  A list with eight data frames:
  \describe{
    \item{\code{ERA5}}{a data frame including the data from reanalysis ERA5.}
    \item{\code{IPSLr1}}{first realization of the IPSL model.}
    \item{\code{IPSLr2}}{second realization of the IPSL model.}
    \item{\code{IPSLr3}}{third realization of the IPSL model.}
    \item{\code{IPSLr4}}{fourth realization of the IPSL model.}
    \item{\code{MIROCr1}}{first realization of the MIROC model.}
    \item{\code{MIROCr2}}{second realization of the MIROC model.}
    \item{\code{MIROCr3}}{third realization of the MIROC model.}
  }
  Each of these data frames includes 4 variables:
  \describe{
    \item{\code{V1}}{a vector with longitude data.}
    \item{\code{V2}}{a vector with latitude data.}
    \item{\code{V3}}{zonal component of surface wind (m/s).}
    \item{\code{V4}}{meridional component of surface wind (m/s).}
  }  

}
\references{
Hersbach H., de Rosnay P., Bell B., Schepers D., Simmons A. et al. (2018). Operational global reanalysis: progress, future directions and synergies with NWP. ECMWF Report. 

Dufresne, J. L., Foujols, M. A., Denvil, S., Caubel, A., Marti, O., Aumont, O., ... & Bony, S. (2013). Climate change projections using the IPSL-CM5 Earth System Model: from CMIP3 to CMIP5. Climate Dynamics, 40(9-10), 2123-2165.
}
\examples{
# Load the data
data(Ensembles)

deg2rad<-function(d){
  return(d*pi/180.)
}

#Get area factors from latitudes
Ws<-sqrt(cos(deg2rad(Ensembles$ERA5$V2)))

#Create a ref/ mod object
#ERA5$V3 is zonal component and ERA5$V4 is meridional component
ref<-data.frame("ERA5",Ensembles$ERA5$V3*Ws,Ensembles$ERA5$V4*Ws)
names(ref)<-c("mod","U","V")

#Create an object for the model runs
modA<-data.frame(mod=rep("MIROCr2",nrow(Ensembles$MIROCr2)),
U=Ensembles$MIROCr2$V3*Ws,V=Ensembles$MIROCr2$V4*Ws)
modB<-data.frame(mod=rep("MIROC",nrow(Ensembles$MIROCr2)),
U=Ensembles$MIROCr2$V3*Ws,V=Ensembles$MIROCr2$V4*Ws)


# In modA every run of the ensemble is considered individually
modA<-rbind(modA,data.frame(mod=rep("MIROCr3",nrow(Ensembles$MIROCr3)),
U=Ensembles$MIROCr3$V3*Ws,V=Ensembles$MIROCr3$V4*Ws))
# In modB all the members of the ensemble are taken as a single 
# model run (too long, so that later on, observations will have to be repeated)
modB<-rbind(modB,data.frame(mod=rep("MIROC",nrow(Ensembles$MIROCr3)),
U=Ensembles$MIROCr3$V3*Ws,V=Ensembles$MIROCr3$V4*Ws))

modA<-rbind(modA,data.frame(mod=rep("MIROCr1",nrow(Ensembles$MIROCr1)),
U=Ensembles$MIROCr1$V3*Ws,V=Ensembles$MIROCr1$V4*Ws))
modB<-rbind(modB,data.frame(mod=rep("MIROC",nrow(Ensembles$MIROCr1)),
U=Ensembles$MIROCr1$V3*Ws,V=Ensembles$MIROCr1$V4*Ws))

modA<-rbind(modA,data.frame(mod=rep("IPSLr4",nrow(Ensembles$IPSLr4)),
U=Ensembles$IPSLr4$V3*Ws,V=Ensembles$IPSLr4$V4*Ws))
modB<-rbind(modB,data.frame(mod=rep("IPSL",nrow(Ensembles$IPSLr4)),
U=Ensembles$IPSLr4$V3*Ws,V=Ensembles$IPSLr4$V4*Ws))

modA<-rbind(modA,data.frame(mod=rep("IPSLr2",nrow(Ensembles$IPSLr2)),
U=Ensembles$IPSLr2$V3*Ws,V=Ensembles$IPSLr2$V4*Ws))
modB<-rbind(modB,data.frame(mod=rep("IPSL",nrow(Ensembles$IPSLr2)),
U=Ensembles$IPSLr2$V3*Ws,V=Ensembles$IPSLr2$V4*Ws))

modA<-rbind(modA,data.frame(mod=rep("IPSLr3",nrow(Ensembles$IPSLr3)),
U=Ensembles$IPSLr3$V3*Ws,V=Ensembles$IPSLr3$V4*Ws))
modB<-rbind(modB,data.frame(mod=rep("IPSL",nrow(Ensembles$IPSLr3)),
U=Ensembles$IPSLr3$V3*Ws,V=Ensembles$IPSLr3$V4*Ws))

modA<-rbind(modA,data.frame(mod=rep("IPSLr1",nrow(Ensembles$IPSLr1)),
U=Ensembles$IPSLr1$V3*Ws,V=Ensembles$IPSLr1$V4*Ws))
modB<-rbind(modB,data.frame(mod=rep("IPSL",nrow(Ensembles$IPSLr1)),
U=Ensembles$IPSLr1$V3*Ws,V=Ensembles$IPSLr1$V4*Ws))

CMIPA<-list(ref,modA)
names(CMIPA)<-c("ref","mod")
CMIPB<-list(ref,modB)
names(CMIPB)<-c("ref","mod")


#Example 1
#Parameters
Uxlim=c(-1,1)
Uylim=c(-0.5,1)
Uxlab<-"U (m/s)"
Uylab<-"V (m/s)"
plotmain<-"CMIP5 samples vs ERA5"
sfactor<-0.2

ref<-CMIPB[["ref"]]
mod<-CMIPB[["mod"]]

sindex<-SailoR.Indices(ref, mod, Ensembles=TRUE)
SailoR.Table(sindex, round_digits = 3)

SailoR.Plot(ref, mod, ColourList=NULL, sfactor, docenter=TRUE,
            Uxlim, Uylim, Uxlab, Uylab, plotmain, plotlegend=TRUE,
            Ensembles=TRUE, plotRMSElegend=TRUE,
            RMSE_legend_Rounding=2, RMSE_legend_units = " m/s",
            referenceName="ERA5")


#Example 2
#Parameters
Uxlim=c(-1,1.25)
Uylim=c(-0.25,1.25)
Uxlab<-"U (m/s)"
Uylab<-"V (m/s)"
plotmain<-"CMIP5 samples vs ERA5"
sfactor<-0.2

ref<-CMIPA[["ref"]]
mod<-CMIPA[["mod"]]


sindex<-SailoR.Indices(ref,mod,Ensembles=TRUE)
SailoR.Table(sindex,round_digits = 3)

SailoR.Plot(ref, mod, ColourList=NULL, sfactor, docenter=FALSE,
            Uxlim, Uylim, Uxlab, Uylab, plotmain, plotlegend=TRUE,
            Ensembles=TRUE, plotRMSElegend=TRUE,
            RMSE_legend_Rounding=2, RMSE_legend_units = " m/s",
            referenceName="ERA5")


}
\keyword{Datasets}
