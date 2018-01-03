\name{ltsk}
\alias{ltsk}
\title{
Ordinary Local Time and Space Kriging
}
\description{
Function implements ordinary time and space kriging for large data sets, with automatic product-sum variogram estimation.
}
\usage{
ltsk(query, obs, th, xcoord = "x", ycoord = "y", tcoord = "t", 
	zcoord = "z", vth = NULL, vlen = NULL, llim = c(3, 3), 
	verbose = T, Large = 2000, future=T, cl = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{query}{
a data frame containing query spatiotemporal locations for which predictions are needed
}
  \item{obs}{
a data frame containing spatiotemporal locations and observed data
}
  \item{th}{
a vector, distance threshold and time lag to define neighbors of a query point
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
	thersholds for local spatiotemporal varigoram (default 75\% of the max lag difference)
}
  \item{vlen}{
	numbers of bins for local spatiotemporal varigram(default,space 15, temporal for each day)
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
	a \pkg{parallel} cluster object (default number of cores in local PC minus one), 0 means single core
}
}
\details{
Function implements automatic variogram estimation (when possible) within a local spatiotemporal neighborhoods, and ordinary krigng based on the produce-sum variogram within that neighborhood. An variogram is estimated for each query point to allow for possible non-stationarity in the data generating field. 

If the number of neighbors exceeds a user-specified upper limit (\code{Large}), neighbors are sub-sampled in a balanced way to reduce the neighborhood size. 

Four variogram models: Gaussian, exponential, spherical and Matern are automatically fit to the empirical space and time variogram in the first lag. The range parameter is estimated from the first distance lag where the empirical variogram exceeds 80\% of the maximum. Weighted least square is then used to estimate the nugget and partial sill parameters. Model with minimal residual sum of squares between the empirical and fitted variogram is chosen as the variogram model.

}
\value{
Kriging mean and standard deviation and quality flags.
\tabular{ll}{
0 \tab valid prediction\cr
1 \tab not enough temporal neighbors \cr
2 \tab not enough spatial neighbors \cr
3 \tab not enough neighbors \cr
4 \tab variogram could not be fit \cr
}

}
\references{
Haas, Timothy C. "Local prediction of a spatio-temporal process with an application to wet sulfate deposition." Journal of the American Statistical Association 90.432 (1995): 1189-1199.

Iaco, S. De & Myers, D. E. & Posa, D., 2001. "Space-time analysis using a general product-sum model," Statistics & Probability Letters, Elsevier, vol. 52(1), pages 21-28, March.

Kumar, N., et al. (2013). "Satellite-based PM concentrations and their application to COPD in Cleveland, OH." Journal of Exposure Science and Environmental Epidemiology 23(6): 637-646.
	
Liang, D. and N. Kumar (2013). "Time-space Kriging to address the spatiotemporal misalignment in the large datasets." Atmospheric Environment 72: 60-69.
}
\author{
Naresh Kumar (NKumar@med.miami.edu)
Dong Liang (dliang@umces.edu)
}

\examples{
## load the data
data(ex)
data(epa_cl)
## apply log transformation
obs[,'pr_pm25'] = log(obs[,'pr_pm25'])
## run kriging
system.time(out <- ltsk(ex2.query[1:2,],obs,c(0.10,10),zcoord='pr_pm25',verbose=FALSE,cl=0))
table(out$flag)
}


