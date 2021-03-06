\name{SailoR-package}
\alias{SailoR-package}
\alias{SailoR}
\docType{package}
\title{
An Extension of the Taylor Diagram to Two-Dimensional Vector Data
}
\description{
A new diagram for the verification of vector variables (wind, current, etc) generated by multiple models against a set of observations is presented in this package. It has been designed as a generalization of the Taylor diagram to two dimensional quantities. It is based on the analysis of the two-dimensional structure of the mean squared error matrix between model and observations. The matrix is divided into the part corresponding to the relative rotation and the bias of the empirical orthogonal functions of the data. The full set of diagnostics produced by the analysis of the errors between model and observational vector datasets comprises the errors in the means, the analysis of the total variance of both datasets, the rotation matrix corresponding to the principal components in observation and model, the angle of rotation of model-derived empirical orthogonal functions respect to the ones from observations, the standard deviation of model and observations, the root mean squared error between both datasets and the squared two-dimensional correlation coefficient. See the output of function \code{UVError()} in this package.
}
\author{
Jon Sáenz, Sheila Carreno-Madinabeitia, Santos J. González-Rojí, 
Ganix Esnaola, Gabriel Ibarra-Berastegi and Alain Ulazia

Maintainer: Santos J. González-Rojí <santosjose.gonzalez@ehu.eus>
}
\keyword{ Package }

\references{
J. Sáenz, S. Carreno-Madinabeitia, G. Esnaola, S.J. González-Rojí, G. Ibarra-Berastegi and A. Ulazia (2020). The Sailor diagram - A new diagram for the verification of two-dimensional vector data from multiple models. Geoscientific Model Development, 13(7), 3221-3240. doi: \url{https://doi.org/10.5194/gmd-13-3221-2020}.
}

\seealso{
 \code{\link{UVError}}
}
\examples{
#--------------------------------------------------------
# Example 1: Figure 5 (left) from Sáenz et al. (2020)
#--------------------------------------------------------

# Load the data
data(Dragonera)

# Parameters
Uxlim=c(-10,10)
Uylim=c(-10,10)
Uxlab<-"U (m/s)"
Uylab<-"V (m/s)"
plotmain<-"Dragonera U10/V10"
sfactor<-1
ref<-Dragonera[["ref"]]
mod<-Dragonera[["mod"]]
    
# Index
SailoR.Indices(ref,mod)

# Plot    
p <- SailoR.Plot(ref,mod,ColourList=NULL,sfactor,docenter=TRUE,
                 Uxlim,Uylim,Uxlab,Uylab,plotmain,plotlegend=TRUE,
                 plotRMSElegend=TRUE,
                 plotscalelegend=TRUE,
                 RMSE_legend_Rounding=1,RMSE_legend_units = " m/s",
                 referenceName="Reference")

#Add segments
segments(-1,-1, x1 = -1, y1 = 1,col="brown",lwd=4)
segments(-1,-1, x1 = 1, y1 = -1,col="brown",lwd=4)
segments(1,-1, x1 = 1, y1 = 1,col="brown",lwd=4)
segments(1,1, x1 = -1, y1 = 1,col="brown",lwd=4)

  }
