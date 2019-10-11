\name{SailoR.Plot}
\alias{SailoR.Plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sailor Diagram
}
\description{
This function generates the Sailor diagram. It calls \code{SailoR.Indices()} internally. 
}
\usage{
SailoR.Plot(ref, mod, ColourList = NULL, sfactor = 1, docenter = FALSE, xlim = NULL, 
ylim = NULL, xlab = NULL, 
ylab = NULL, plotmain = NULL, plotlegend = TRUE, 
plotRMSElegend = TRUE, Ensembles = FALSE, RMSE_legend_units = "", 
RMSE_legend_Rounding = 2, referenceName = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ref}{
a data frame with the reference observations. It should have 3 columns: name, zonal component and meridional component.
}
  \item{mod}{
a data frame with the model data. It should have 3 columns: name, zonal component and meridional component. More than one model can be included by using different names. 
}
  \item{ColourList}{
a vector with all the colors you want to use for plotting the models. Grey will always be used for reference data. If \code{NULL}, the default colours are used: c("darkgrey", "red", "blue", "darkolivegreen3", "orange", "seagreen3", "gold", "purple", "pink").}
  \item{sfactor}{
a value with the scale factor to be applied to the ellipses of each model.}
  \item{docenter}{
if all the ellipses should be plotted over the reference one, this argument should be \code{TRUE}.}
  \item{xlim}{
X axis limit for the plot.
}
  \item{ylim}{
Y axis limit for the plot.
}
  \item{xlab}{
a title for the X axis.
}
  \item{ylab}{
a title for the Y axis.
}
  \item{plotmain}{
an overall title for the plot.
}
  \item{plotlegend}{
if a default legend should be added to the plot, this argument should be \code{TRUE}.
}
  \item{plotRMSElegend}{
if a legend with a summary of the RMSE values should be added to the plot, this argument should be \code{TRUE}.
}
  \item{Ensembles}{
if \code{Ensembles} is set to \code{TRUE}, M different models included in the data frame \code{mod} will be taken as a single model and the reference model is repeated M times for the calculation of the Sailor diagram.
}
  \item{RMSE_legend_units}{
units to be written in the RMSE legend.
}
  \item{RMSE_legend_Rounding}{
number of decimals to be kept in the RMSE legend.
}
  \item{referenceName}{
external name provided for the reference dataset for the legend.
}
}
\value{
The result is a plot object.
}
\seealso{
 \code{\link{SailoR.Indices}}
}
\examples{
# Load the Vertically Integrated Water Vapor Transport Data included in the package
data(WRF)

#Parameters
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

#Index
sIWRF=SailoR.Indices(ref,mod)

#Index table
sIWRF.table<-SailoR.Table(sIWRF)

#Save plot
SailoR.Plot(ref,mod,ColourList=NULL,sfactor,docenter=TRUE,
            Uxlim, Uylim, Uxlab, Uylab, plotmain, plotlegend=TRUE,
            plotRMSElegend=TRUE, RMSE_legend_Rounding=0,
            RMSE_legend_units=" kg/m/s", referenceName="Sounding")
  }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Functions}% use one of  RShowDoc("KEYWORDS")