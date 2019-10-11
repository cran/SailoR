\name{Dragonera}
\alias{Dragonera}
\docType{data}
\title{
Wind data
}
\description{
Contains the information of both zonal and meridional wind components measured in the Mediterranean sea by the in-situ buoy Dragonera, operated by the Spanish State Ports Authority - Puertos del Estado. Data from two WRF downscaling experiments (including the 3DVAR data assimilation every 3 hours - D, and not including it - N), ERA-Interim reanalysis and remote sensing data from the second version of the Cross-Calibrated Multi-platform (CCMPv2) product are also included. Only the closest spatial grid points to the buoy location (39.56 ºN, 2.10 ºE) in the gridded products are considered. The data cover the period 2009-2014 (Ulazia et al., 2017). 
}
\usage{data("Dragonera")}
\format{
   A list with two data frames:
  \describe{
    \item{\code{ref}}{a data frame including the buoy data.}
    \item{\code{mod}}{a data frame including the data from gridded datasets.}
  } 
  Each of these data frames includes 3 variables:
  \describe{
    \item{\code{mod}}{a factor defining if the data belong to the CCMPv2 product (CCMP_SAT), to the reanalysis ERA-Interim (rad) or to any of the WRF downscaling experiments (WRF_D or WRF_N).}
    \item{\code{U}}{zonal component of wind (m/s).}
    \item{\code{V}}{meridional component of wind (m/s).}
  }  
}

\source{
Both WRF downscaling experiments were created by the authors in previous studies. The in-situ data can be retrieved from \url{http://www.puertos.es/es-es/oceanografia/Paginas/portus.aspx}. The CCMPv2 data can be retrieved from \url{http://www.remss.com/measurements/ccmp/}
}

\references{
A. Ulazia, J. Sáenz, G. Ibarra-Berastegi, S. J. González-Rojí and S. Carreno-Madinabeitia (2017) Using 3DVAR data assimilation to measure offshore wind energy potential at different turbine heights in the West Mediterranean. Applied Energy, 208:1232-1245, doi: \url{https://doi.org/10.1016/j.apenergy.2017.09.030}
}

\examples{
# Load the data
data(Dragonera)

#Example 1
#Parameters
Uxlim=c(-10,10)
Uylim=c(-10,10)
Uxlab<-"U (m/s)"
Uylab<-"V (m/s)"
plotmain<-"Dragonera U10/V10"
sfactor<-1
ref<-Dragonera[["ref"]]
mod<-Dragonera[["mod"]]
    
#Index
SailoR.Indices(ref,mod)
    
SailoR.Plot(ref, mod, ColourList=NULL, sfactor, docenter=FALSE,
            Uxlim, Uylim, Uxlab, Uylab, plotmain, plotlegend=TRUE,
            plotRMSElegend=TRUE, RMSE_legend_Rounding=1, 
            RMSE_legend_units = " m/s", referenceName="Reference")
  
#Example 2
#Parameters
Uxlim=c(-1,1)
Uylim=c(-1,1)
Uxlab<-"U (m/s)"
Uylab<-"V (m/s)"
plotmain<-"Dragonera U10/V10"
sfactor<-0.025
ref<-Dragonera[["ref"]]
mod<-Dragonera[["mod"]]
  
#Index
SailoR.Indices(ref, mod)
  
plotmain<-"Dragonera U10/V10 (scaled)"
SailoR.Plot(ref, mod, ColourList=NULL, sfactor, docenter=FALSE,
            Uxlim, Uylim, Uxlab, Uylab, plotmain, plotlegend=TRUE,
            plotRMSElegend=TRUE, RMSE_legend_Rounding=1,
            RMSE_legend_units=" m/s", referenceName="Reference")

}
\keyword{Datasets}
