\name{Reanalysis}
\alias{Reanalysis}
\docType{data}
\title{
Surface Wind from different Reanalyses
}
\description{
Contains the information of both components of monthly averaged surface wind in January over the Northern Hemisphere as included in different state-of-the-art reanalyses during the period 2011-2018. ERA5 reanalysis was defined as the reference, and the remaining gridded products were included in the model matrix. Particularly, the original NCEP/NCAR first generation reanalysis, MERRA2, CFSv2 and ERA-Interim. All of them have been bilinearly interpolated to NCEP/NCAR original grid.
}
\usage{data("Reanalysis")}
\format{
   A list with two matrices:
  \describe{
    \item{\code{ref}}{a matrix including the data from reanalysis ERA5.}
    \item{\code{mod}}{a matrix including the data from other reanalyses.}
  } 
  Each of these data frames includes 5 variables:
  \describe{
    \item{\code{V1}}{a vector with longitude data.}
    \item{\code{V2}}{a vector with latitude data.}
    \item{\code{V3}}{a factor defining if the data belong to the reanalysis CFSv2 (cfsv2), ERA-Interim (ei), MERRA2 (merra2) or NCEP/NCAR reanalysis (nnra).}
    \item{\code{V4}}{zonal component of surface wind (m/s).}
    \item{\code{V5}}{meridional component of surface wind (m/s).}
      }  
}
\examples{
# Load the data
data(Reanalysis)

#Parameters
Uxlim=c(-1,1)
Uylim=c(-1,1)
Uxlab<-"U (m/s)"
Uylab<-"V (m/s)"
plotmain<-"Reanalyses"
sfactor<-0.2

#Create ref, mod objects

weightWithLat=TRUE

if(weightWithLat){
    dlats=as.numeric(Reanalysis[["ref"]][,2])
    rlats=pi*dlats/180.
    lweights=sqrt(cos(rlats))

    ref<-data.frame(Reanalysis[["ref"]][,3],
                      lweights*as.numeric(Reanalysis[["ref"]][,4]),
                      lweights*as.numeric(Reanalysis[["ref"]][,5]))
}else{

    ref<-data.frame(Reanalysis[["ref"]][,3],
                      as.numeric(Reanalysis[["ref"]][,4]),
                      as.numeric(Reanalysis[["ref"]][,5]))

}
names(ref)<-c("mod","U","V")

mod<-data.frame(Reanalysis[["mod"]][,3],
                  as.numeric(Reanalysis[["mod"]][,4]),
                  as.numeric(Reanalysis[["mod"]][,5]))
names(mod)<-c("mod","U","V")
  
#Sailor

SailoR.Plot(ref, mod, ColourList=NULL, sfactor, docenter=TRUE,
              Uxlim, Uylim, Uxlab, Uylab, plotmain, plotlegend=TRUE,
              plotRMSElegend=TRUE, RMSE_legend_Rounding=1,
              RMSE_legend_units = " m/s", referenceName="ERA5")


SailoR.Plot(ref, mod, ColourList=NULL, sfactor, docenter=TRUE,
              Uxlim, Uylim, Uxlab, Uylab, plotmain, plotlegend=TRUE,
              plotRMSElegend=TRUE)
  
}
\keyword{Datasets}
