\name{SailoR.Table}
\alias{SailoR.Table}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary of indices calculated for the Sailor diagrams
}
\description{
This function reads the output from \code{SailoR.Indices()} and creates a summary data frame of the indices calculated for the Sailor diagram. This can easily be converted to a LaTeX or other formats with external packages (\code{'xtable'}, for instance).
}
\usage{
SailoR.Table(UV, round_digits = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{UV}{
a list containing the output obtained directly from \code{SailoR.Indices()} function.
}
  \item{round_digits}{
number of decimals to be kept in the summary table.
}
}
\value{
The result is a data frame with the information from the \code{SailoR.Indices()} function in tabular format. Particularly, the indices stated in the table are the following ones:

\item{modelName}{the names of the models datasets.}
\item{sdUx}{standard deviation of zonal component (U).}
\item{sdUy}{standard deviation of meridional component (U).}
\item{sdVx}{standard deviation of zonal component (V).}
\item{sdVy}{standard deviation of meridional component (V).}
\item{Sigmax}{a scalar representing the semi-major axis of the ellipse (standard deviation of PC1).}
\item{Sigmay}{a scalar representing the semi-minor axis of the ellipse (standard deviation of PC2).}
\item{thetau}{angle of EOF1 (U) with zonal axis.}
\item{thetav}{angle of EOF1 (V) with meridional axis.}
\item{thetavu}{rotation angle from U EOFs to V EOFs.}
\item{R2vec}{vector correlation squared in 2D as defined by Breaker, Gemmill and Crossby (1994).}
\item{biasMag}{magnitude of the bias vector.}
\item{RMSE}{root mean square error between U and V.}
\item{Eccentricity}{eccentricity of the ellipses from the reference and models dataset.}
\item{congruenceEOF1}{congruence coefficients (absolute value) for EOF1.}

}
\examples{
#--------------------------------------------------------
# Example 1
#--------------------------------------------------------
# Load the ocean current data 
data("Current")
  
# Calculate the statistics
ref<-Current[["ref"]]
mod<-Current[["mod"]]
UV <- SailoR.Indices(ref,mod)

# Obtain the summary of the results
SailoR.Table(UV)

#--------------------------------------------------------
# Example 2: Table 2 from Sáenz et al. (2020)
#--------------------------------------------------------
# Load Synthetic data 
data(Synthetic)

ref<-Synthetic[["ref"]]
mod<-Synthetic[["mod"]]

# Calculate the statistics
UVtable=SailoR.Table(UV,round_digits = 4)
TotVar<-c((UVtable$sdUx^2+UVtable$sdUy^2)[1],
          (UVtable$sdVx^2+UVtable$sdVy^2)[-1])
TotSigma<-UVtable$Sigmax^2+UVtable$Sigmay^2
UVtable<-cbind(UVtable[,1],TotVar,TotSigma,UVtable[,8:15])

UVtable<-round(UVtable[,-1],2)

#If you want to print it as a LaTeX table:
#library(xtable)
#oInfo<-xtable(UVtable)
#print(oInfo,type="latex",file="table.tex")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Functions}% __ONLY ONE__ keyword per line
