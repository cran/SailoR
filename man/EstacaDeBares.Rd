\name{EstacaDeBares}
\alias{EstacaDeBares}
\docType{data}
\title{
Wave Energy Flux data
}
\description{
Contains the information of both components of wave energy flux measured by the buoy placed in Estaca de Bares and operated by the Spanish State Ports Authority - Puertos del Estado. Different statistical models were applied to this observed data: analogs, random forests, the combination of both of them, WAM model and the persistence (Ibarra-Berastegi et al., 2015 and 2016).   
}
\usage{data("EstacaDeBares")}
\format{
 A list with two data frames:
  \describe{
    \item{\code{ref}}{a data frame including the buoy data.}
    \item{\code{mod}}{a data frame including the data from the statistical models.}
  } 
  Each of these data frames includes 4 variables:
  \describe{
    \item{\code{date}}{date information.}
    \item{\code{mod}}{a factor defining if the data belong to the analogs method (analo), to random forests (rf), both techniques combined (analorf), WAM model (wam) or persistence (pers).}
    \item{\code{U}}{zonal component of wave energy flux (kW/m).}
    \item{\code{V}}{meridional component of wave energy flux (kW/m).}
  }    
}
\source{
The results from the statistical methods were obtained by the authors in previous studies. The in-situ data can be retrieved from \url{http://www.puertos.es/es-es/oceanografia/Paginas/portus.aspx}.
}

\references{
G. Ibarra-Berastegi, J. Sáenz, G. Esnaola, A. Ezcurra and A. Ulazia (2015) Short-term forecasting of the wave energy flux: Analogues, random forests, and physics-based models. Ocean Engineering, 104:530-539, doi: \url{https://doi.org/10.1016/j.oceaneng.2015.05.038}

G. Ibarra-Berastegi, J. Sáenz, G. Esnaola, A. Ezcurra, A. Ulazia, N. Rojo and G. Gallastegui (2016) Wave Energy Forecasting at Three Coastal Buoys in the Bay of Biscay IEEE JOURNAL OF OCEANIC ENGINEERING, 41:923-929, doi: \url{https://doi.org/10.1109/JOE.2016.2529400}
}

\examples{
# Load the data
data(EstacaDeBares)

#Parameters
Uxlim=c(0,25)
Uylim=c(0,25)
Uxlab<-" "
Uylab<-" "
plotmain<-"Wave Energy Flux. Estaca De Bares."
sfactor<-0.25
plotRMSE<-TRUE
ref<-EstacaDeBares[["ref"]][,c("mod","U","V")]
mod<-EstacaDeBares[["mod"]][,c("mod","U","V")]
  
#Index
SailoR.Indices(ref,mod)
  

SailoR.Plot(ref, mod, ColourList=NULL, sfactor, docenter=TRUE,
            Uxlim, Uylim, Uxlab, Uylab, plotmain, plotlegend=TRUE,
            plotRMSElegend=TRUE, RMSE_legend_Rounding=0, 
            RMSE_legend_units=" W/m", referenceName="Buoy")

}
\keyword{Datasets}
