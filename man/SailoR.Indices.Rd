\name{SailoR.Indices}
\alias{SailoR.Indices}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculation of indices for Sailor Diagram
}
\description{
This function returns a data frame with all the statistics used by the Sailor diagram for the set of different models in data frame \code{mod} with the reference in data frame \code{ref}.
}
\usage{
SailoR.Indices(ref, mod, Ensembles = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ref}{
A data frame with the reference observations. It should have 3 columns: name, zonal component and meridional component.
}
  \item{mod}{
A data frame with the model data. It should have 3 columns: name, zonal component and meridional component. More than one model can be included under different names. Additionally, \code{nrows(mod)} must be an integer multiple of \code{nrows(ref)}. 
}
  \item{Ensembles}{
If \code{Ensembles} is set to \code{TRUE}, M different models included in the data frame \code{mod} will be taken as a single model and the reference model is repeated M times for the calculation of the Sailor diagram. 
}
}
\value{
The result is a data frame with a summary of the statistics obtained:

\item{meanU}{mean of the U dataset.}
\item{meanV}{mean of the V dataset.}
\item{RMSE}{root mean square error between U and V.}
\item{sdUx}{standard deviation of zonal component (U).}
\item{sdUy}{standard deviation of meridional component (U).}
\item{sdVx}{standard deviation of zonal component (V).}
\item{sdVy}{standard deviation of meridional component (V).}
\item{Eu}{EOF matrix for the U dataset.}
\item{Ev}{EOF matrix for the V dataset.}
\item{Sigmau}{matrix with standard deviations of U (2x2 matrix).}
\item{Sigmav}{matrix with standard deviations of V (2x2 matrix).}
\item{TotVarU}{total variance of the U dataset.}
\item{TotVarV}{total variance of the V dataset.}
\item{thetau}{angle of EOF1 (U) with zonal axis.}
\item{thetav}{angle of EOF1 (V) with meridional axis.}
\item{thetavu}{rotation angle from U EOFs to V EOFs.}
\item{R2vec}{vector correlation squared in 2D as defined by Breaker, Gemmill and Crossby (1994).}
\item{Rvu}{rotation matrix to express Ev as a rotation from Eu.}
\item{EccentricityU}{eccentricity of the ellipses from the reference dataset.}
\item{EccentricityV}{eccentricity of the ellipses from the model dataset.}
\item{congruenceEOF}{congruence coefficients (absolute value) for EOFs.}

}

\references{
L. C. Breaker, W. H. Gemmill and D. S. Crosby (1994). The application of a technique for vector correlation to problems in meteorology and oceanography. Journal of Applied Meteorology, 33(11), 1354-1365.
}

\examples{
#--------------------------------------------------------
# Example 1
#--------------------------------------------------------

# Load the ocean current data 
data("Current")
  
# Calculate the indices
ref<-Current[["ref"]]
mod<-Current[["mod"]]

SailoR.Indices(ref, mod)
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Functions}% use one of  RShowDoc("KEYWORDS")
