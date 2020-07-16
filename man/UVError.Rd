\name{UVError}
\alias{UVError}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculation of the elements of the 2D MSE terms
}
\description{
This function calculates the different terms in the 2D MSE equation, the rotation of EOFs and so on for a unique given model. This function is internally used by \code{SailoR.Indices()} and external users are not expected to use it.
}
\usage{
UVError(U, V, Ensembles = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{U}{
a vector including the zonal and meridional components of the reference dataset.}
  \item{V}{
a vector including the zonal and meridional components of the model being compared.}
  \item{Ensembles}{
if \code{Ensembles} is set to \code{TRUE}, M different models included in the data frame \code{mod} will be taken as a single model and the reference model is repeated M times for the calculation of the Sailor diagram. 
}
}
\value{
The result is a list with a summary of the statistics obtained:

\item{meanU}{mean of the U dataset.}
\item{meanV}{mean of the V dataset.}
\item{TotVarU}{total variance of the U dataset.}
\item{TotVarV}{total variance of the V dataset.}
\item{Eu}{EOF matrix for the U dataset.}
\item{Ev}{EOF matrix for the V dataset.}
\item{Rvu}{rotation matrix to express Ev as a rotation from Eu.}
\item{Sigmau}{matrix with standard deviations of U (2x2 matrix).}
\item{Sigmav}{matrix with standard deviations of V (2x2 matrix).}
\item{sdUx}{standard deviation of zonal component (U).}
\item{sdUy}{standard deviation of meridional component (U).}
\item{sdVx}{standard deviation of zonal component (V).}
\item{sdVy}{standard deviation of meridional component (V).}
\item{thetau}{angle of EOF1 (U) with zonal axis.}
\item{thetav}{angle of EOF1 (V) with meridional axis.}
\item{thetavu}{rotation angle from U EOFs to V EOFs.}
\item{RMSE}{root mean square error between U and V.}
\item{R2vec}{vector correlation squared in 2D as defined by Breaker, Gemmill and Crossby (1994).}
\item{EccentricityU}{eccentricity of the ellipses from the reference dataset.}
\item{EccentricityV}{eccentricity of the ellipses from the model dataset.}
\item{congruenceEOF}{congruence coefficients (absolute value) for EOFs.}

}

\references{
L. C. Breaker, W. H. Gemmill, and D. S. Crosby (1994). The application of a technique for vector correlation to problems in meteorology and oceanography. Journal of Applied Meteorology, 33(11), 1354-1365.
}

\seealso{
 \code{\link{SailoR.Indices}}
}

\examples{
# No examples are given as external users are not expected to use it
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Functions}% use one of  RShowDoc("KEYWORDS")

