\name{cltsk}
\alias{cltsk}
\title{
Function calls \code{ltsk} using cumulatively expanding time space thresholds. This function is useful when predictions are needed using data points at different spatiotemporal intervals. 
For example, if predictions are needed at a given location for the past 30 days at an interval of 3 days. 
Instead of using \code{ltsk} 10 times, \code{cltsk} can compute all 10 values simultaneously.
}
\description{
Function calls \code{ltsk} using cumulatively expanding time space thresholds.}
\usage{
cltsk(query, obs, th, nbins, xcoord = "x", ycoord = "y", tcoord = "t", 
	zcoord = "z", vth = NULL, vlen = NULL, llim = c(3, 3), 
	verbose = T, Large = 2000, future=T,cl = NULL)
}
\arguments{
  \item{query}{
data frame containing query point (X,Y,T i.e. XY coordinates and time) where predictions are needed
}
  \item{obs}{
data frame containing sample data with XY coordinates, time and observed (measured) values
}
  \item{th}{
a priori chosen distance and time thresholds for neighbor search
}
  \item{nbins}{
  a vector, number of distance and time bins for cumulative neighbor search and kriging.
  }
  \item{xcoord}{
a character constant, the field name for x coordinate in both \code{query} and \code{obs}
}
  \item{ycoord}{
a character constant, the field name for y coordinate in both \code{query} and \code{obs}
}
  \item{tcoord}{
a character constant, the field name for time coordinate in both \code{query} and \code{obs}
}

  \item{zcoord}{
a character constant, the field name for data in \code{obs}
}
  \item{vth}{
	thresholds for local spatiotemporal variogram (default 75\% of the max lag difference)
}
  \item{vlen}{
	numbers of bins for local spatiotemporal variogram(default, space 15, temporal for each day)
}
  \item{llim}{
lower limits for number of regions and intervals with observed data to calculate Kriging (default 3 spatial regions, 3 temporal intervals)
}
  \item{verbose}{
logical, whether print details information
}
  \item{Large}{
a numeric constant, upper limit of neighbor points, beyond which subsampling is performance
}
	\item{future}{
logical,	whether including observed points in future relative to query points.
}
  \item{cl}{
	a \pkg{parallel} cluster object (default single core)
}

}
\details{
Function performs automatic variogram estimation for each query location using the observed data within \code{th} thresholds. The estimated variogram is used for ordinary kriging, but using data in expanding local neighborhoods for ordinary kriging. 
For example, if predictions are needed at a given location for the past 30 days at an interval of 3 days,data within 3 days are used first, followed by 6 days and so on until data within 30 days. The same applies for distance thresholds.
}
\value{
\enumerate{
\item \verb{krig} Kriging estimates at each space and time neighborhood
\item \verb{legend} The legend for space and time neighborhood
}
}
\references{
Iaco, S. De & Myers, D. E. & Posa, D., 2001. "Space-time analysis using a general product-sum model," Statistics & Probability Letters, Elsevier, vol. 52(1), pages 21-28, March.

Kumar, N., et al. (2013). "Satellite-based PM concentrations and their application to COPD in Cleveland, OH." Journal of Exposure Science and Environmental Epidemiology 23(6): 637-646.
	
Liang, D. and N. Kumar (2013). "Time-space Kriging to address the spatiotemporal misalignment in the large datasets." Atmospheric Environment 72: 60-69.
}
\author{
Naresh Kumar (nkumar@med.miami.edu)
Dong Liang (dliang@umces.edu)
}

\examples{
## load the data
data(ex)
data(epa_cl)
## apply log transformation
obs[,'pr_pm25'] = log(obs[,'pr_pm25'])
## run kriging
system.time(out <- cltsk(ex2.query[1:2,],obs,c(0.10,10),
  zcoord='pr_pm25',nbins=c(4,5),verbose=FALSE))
table(out$flag)
}


