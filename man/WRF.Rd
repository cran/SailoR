\name{WRF}
\alias{WRF}
\docType{data}
\title{
Vertically Integrated Water Vapor Transport Data
}
\description{
Contains the information of both components of the vertically integrated water vapor transports calculated over A Coruña (43.6 ºN, -8.41 ºE) in the Iberian Peninsula. This quantity is obtained after integrating the product of moisture and both components of wind across the whole vertical dimension of the atmosphere. To do so, moisture and wind components measured by radiosondes were retrieved from the server of the University of Wyoming, and used for its calculation. Data from two WRF downscaling experiments (including the 3DVAR data assimilation every 3 hours - D experiment, and not including it - N experiment) and ERA-Interim data were also used for its calculation. In these cases, only the closest spatial grid points to the station in the gridded products are considered. 12-hourly data from period 2010-2014 are included with missing soundings removed (González-Rojí et al., 2018 and 2019). 
}
\usage{data("WRF")}
\format{
  A list with two data frames:
  \describe{
    \item{\code{ref}}{a data frame including the vertically integrated water vapor transports calculated by means of radiosonde data.}
    \item{\code{mod}}{a data frame including the data from gridded datasets.}
  } 
  Each of these data frames includes 4 variables:
  \describe{
    \item{\code{date}}{date information.}
    \item{\code{mod}}{a factor defining if the values were calculated by means of reanalysis ERA-Interim data (ERAI) or data from WRF downscaling experiments (WrfD or WrfN).}
    \item{\code{U}}{zonal component of vertically integrated water vapor transport (kg/m/s).}
    \item{\code{V}}{meridional component of vertically integrated water vapor transport (kg/m/s).}
  }  
}
\source{
Both downscaling experiments were created by the authors in previous studies. ERA-Interim were retrieved from the MARS repository at ECMWF.  
}
\references{
S. J. González-Rojí, J. Sáenz, G. Ibarra-Berastegi and J. Díaz-Argandoña (2018) Moisture balance over the Iberian Peninsula according to a regional climate model. The impact of 3DVAR data assimilation. Journal of Geophysical Research-Atmospheres, 123:708-729, doi: \url{https://doi.org/10.1002/2017JD027511}

S. J. González-Rojí, R. L. Wilby, J. Sáenz and G. Ibarra-Berastegi (2019) Harmonized evaluation of daily precipitation downscaled using SDSM and WRF+WRFDA models over the Iberian Peninsula Climate Dynamics, 53:1413-1433, doi: \url{https://doi.org/10.1007/s00382-019-04673-9}
}
\examples{
#--------------------------------------------------------
# Example 1: Figure 6 (right) from Sáenz et al. (2020)
#--------------------------------------------------------

# Load the data
data("WRF")

# Parameters
Uxlim=c(60,110)
Uylim=c(0,50)
Uxlab<-"Qx (kg/m/s)"
Uylab<-"Qy (kg/m/s)"
plotmain<-"Water vapour transport"
sfactor<-0.1
ref<-WRF[["ref"]][,c("mod","U","V")]
mod<-WRF[["mod"]][,c("mod","U","V")]

isBad=((is.na(ref$U)) |(is.na(ref$V)))
isOK=(!isBad)
ref<-ref[isOK,]

nmod1<-which(mod$mod=="WrfN")
mod1<-mod[nmod1,][isOK,]
nmod2<-which(mod$mod=="WrfD")
mod2<-mod[nmod2,][isOK,]
nmod3<-which(mod$mod=="ERAI")
mod3<-mod[nmod3,][isOK,]

mod<-rbind(mod1,mod2,mod3)

# Index
sIWRF=SailoR.Indices(ref,mod)

# Index table
sIWRF.table<-SailoR.Table(sIWRF)

# plot
      
SailoR.Plot(ref, mod, ColourList=NULL, sfactor, docenter=TRUE,
            Uxlim, Uylim, Uxlab, Uylab, plotmain, plotlegend=TRUE,
            plotRMSElegend=TRUE, plotscalelegend=TRUE,
            RMSE_legend_Rounding=0, RMSE_legend_units = " kg/m/s",
            referenceName="Sounding")
}
\keyword{Datasets}
