\name{Synthetic}
\alias{Synthetic}
\docType{data}
\title{
Ocean Current Data
}
\description{
This dataset contains a one-year long time-series of hourly zonal and meridional wind components in the extratropical Pacific from ERA5. From this dataset (Reference), different synthetic datasets with individual sources of error artificially added to them have been produced. In the first case (MOD1), a constant bias has been added. Next, the error is produced (MOD2) by rotating the zonal and meridional components by 30 degrees counterclockwise (thus also inducing a bias due to the rotation of the mean vector). An unphysical source of error is added in MOD3 by randomly resampling the dataset in order to break the original correlation of the vectors while keeping the bias and EOFs at their original values. Finally, the original data has been multiply by a constant in order to change the variance of the data (MOD4), although the scaling produces a change in the bias, too. The main objective of adding this synthetic dataset is to show the response of the Sailor diagram to different kinds of errors.
}
\usage{data("Synthetic")}
\format{
  A list with two data frames:
  \describe{
    \item{\code{mod}}{a data frame including the wind data from ERA5.}
    \item{\code{ref}}{a data frame including the synthetic datasets.}
  } 
  Each of these data frames includes 3 variables:
  \describe{
    \item{\code{mod}}{a factor defining if the data belong to the ERA5 wind data (REF), or to the synthetic datasets MOD1, MOD2, MOD3 or MOD4.}
    \item{\code{U}}{zonal component of wind speed (m/s).}
    \item{\code{V}}{meridional component of wind speed (m/s).}
  }  
}

\examples{
# Load the data
data(Synthetic)

Uxlim=c(0,15)
Uylim=c(-25,10)
Uxlab<-"Ux (m/s)"
Uylab<-"Uy (m/s)"
plotmain<-"Reference and synthetic models"
sfactor<-1
ref<-Synthetic[["ref"]]
mod<-Synthetic[["mod"]]

#--------------------------------------------------------
# Example 1: Figure 4 (left) from Sáenz et al. (2020)
#--------------------------------------------------------

SailoR.Plot(ref, mod, ColourList=NULL, sfactor, docenter=FALSE,
            Uxlim, Uylim, Uxlab, Uylab, plotmain, plotlegend=TRUE,
            Ensembles=TRUE, plotRMSElegend=TRUE, 
            plotscalelegend=TRUE, RMSE_legend_Rounding=2,
            RMSE_legend_units = " m/s")

#--------------------------------------------------------
# Example 2: Figure 4 (right) from Sáenz et al. (2020)
#--------------------------------------------------------

SailoR.Plot(ref, mod, ColourList=NULL, sfactor, docenter=TRUE,
            Uxlim, Uylim, Uxlab, Uylab, plotmain, plotlegend=TRUE,
            Ensembles=TRUE, plotRMSElegend=TRUE,
            plotscalelegend=TRUE, RMSE_legend_Rounding=2,
            RMSE_legend_units = " m/s")
}

\keyword{Datasets}
