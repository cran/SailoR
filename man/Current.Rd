\name{Current}
\alias{Current}
\docType{data}
\title{
Ocean Current Data
}
\description{
Contains the information of the surface ocean currents measured in the Bay of Biscay by the Donostia buoy and the HF-Radar from the Basque Meteorological Agency EUSKALMET. Additionally, it also contains model data from the high resolution global analysis and forecasting system PSY4V3R1 of the Copernicus Marine Environment Monitoring Service (CMEMS). Only the closest spatial grid points to the buoy location (43.6 ºN and 2.0 ºW) in radar and model grids are considered, and only the parts of years 2017 and 2018 which are present in the three sources (8319 hourly observations) are included. 
}
\usage{data("Current")}
\format{
  A list with two data frames:
  \describe{
    \item{\code{mod}}{a data frame including the data from both gridded datasets.}
    \item{\code{ref}}{a data frame including the buoy data.}
  } 
    Each of these data frames includes 3 variables:
  \describe{
    \item{\code{mod}}{a factor defining if the data belong to the observed data (ref), to the radar (rad) or to the ocean modelling product (mod).}
    \item{\code{U}}{zonal component of ocean current (m/s).}
    \item{\code{V}}{meridional component of ocean current (m/s).}
  }  
}

\source{
The in-situ data belong to the Basque Meteorological Agency EUSKALMET and were downloaded from \url{https://www.euskoos.eus}. The Copernicus Marine Environment Monitoring Service (CMEMS) data can be retrieved from their data portal.
}

\examples{
# Load the data
data(Current)

#Parameters
Uxlim=c(-0.5,0.5)
Uylim=c(-0.5,0.5)
Uxlab<-"m/s"
Uylab<-"m/s"
plotmain<-"Surface current"
sfactor<-4
ref<-Current[["ref"]]
mod<-Current[["mod"]]
  
#Index
SailoR.Indices(ref,mod)
  
SailoR.Plot(ref, mod, ColourList=NULL, sfactor, docenter=TRUE,
            Uxlim, Uylim, Uxlab, Uylab, plotmain, plotlegend=TRUE,
            plotRMSElegend=TRUE, RMSE_legend_Rounding=2,
            RMSE_legend_units=" m/s", referenceName="Buoy")
}
\keyword{Datasets}
