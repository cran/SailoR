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
The result is a data frame with the information from the \code{SailoR.Indices()} function in table format.
}
\examples{
# Load the ocean current data 
data("Current")
  
# Calculate the Statistics
ref<-Current[["ref"]]
mod<-Current[["mod"]]
UV <- SailoR.Indices(ref,mod)

# Obtain the summary of the results
SailoR.Table(UV)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Functions}% __ONLY ONE__ keyword per line