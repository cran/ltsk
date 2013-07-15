\name{ltsk}
\alias{ltsk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Local time and space kriging
}
\description{
This function implements spatiotemporal neighbor search, automatic variogram estimation using the product-sum model and ordinary Kriging at multiple query locations using the observed data.
}
\usage{
ltsk(query, obs, th, xcoord = "x", ycoord = "y", tcoord = "t", 
	zcoord = "z", vth = NULL, vlen = NULL, llim = c(3, 3), 
	verbose = T, Large = 2000, nproc = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{query}{
data frame containing query spatiotemporal locations
}
  \item{obs}{
data frame containing spatiotemporal locations and observed data
}
  \item{th}{
a priori chosen distance and time thresholds for neighbor search
}
  \item{xcoord}{
field name for x coordinate in both query and obs
}
  \item{ycoord}{
field name for y coordinate in both query and obs
}
  \item{tcoord}{
	field name for time coordinate in both query and obs
}

  \item{zcoord}{
	field name for data in obs
}
  \item{vth}{
	thersholds for local spatiotemporal varigoram (default 75\% max lag difference)
}
  \item{vlen}{
	numbers of bins for local spatiotemporal varigram(default,space 15, temporal for each day)
}
  \item{llim}{
lower limits for number of data points to calculate Kriging (default 3 spatial, 3 temporal neighbors)
}
  \item{verbose}{
boolean whether print details information
}
  \item{Large}{
upper limit of neighbor points, beyond which subsampling was done
}
  \item{nproc}{
	number of cores to use (defautl single core)
}
}
\details{
This function generalizes the local Kriging method to large spatiotemporal data. For each query point, a local 
spatiotemporal neighborhood of observed data was identified. Local spatiotemporal variogram was estiamted for ordinary Kriing. 

Details on the development of spatiotemporal neighborhood through sub-sampling can be found in the local spatial Kriging documents. The product-sum approach wass used to model the local variogram. Details on estimating the variogram model within each neighborhood can be found in Iaco et. al 2001.
}
\value{
Kriging mean and standard error and quality flags.
\tabular{ll}{
0 \tab valid prediction\cr
1 \tab not enough temporal neighbours (less than 3 by default)\cr
2 \tab not enough spatial neighbours (less than 3 by default) \cr
4 \tab not enough neighbours (less than 5 observed data) \cr
}

}
\references{
Iaco, S. De & Myers, D. E. & Posa, D., 2001. "Space-time analysis using a general product-sum model," Statistics & Probability Letters, Elsevier, vol. 52(1), pages 21-28, March.

Kumar,N., Liang, D., Chu, A., Abrams, T., forthcoming. "Satellite based PM concentrations and its application to COPD in Cleveland", Journal of Exposure Science and Environmental Epidemiology.
}
\author{
Dong Liang (dong-liang@uiowa.edu)
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## load the data
data(ex)
data(epa_cl)
## apply log transformation
obs[,'pr_pm25'] = log(obs[,'pr_pm25'])
## run kriging
system.time(out <- ltsk(ex2.query[1:2,],obs,c(0.10,10),zcoord='pr_pm25'))
## simple parallel kriging
system.time(out <- ltsk(ex2.query[1:2,],obs,c(0.10,10),zcoord='pr_pm25',nproc=2))
## output
table(out$flag)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

