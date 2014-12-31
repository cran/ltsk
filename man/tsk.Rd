\name{tsk}
\alias{tsk}
\title{
Ordinary Global Time and Space Kriging
}
\description{
Function for ordinary kriging in time and space based on the product-sum variogram model, kriging in a local neighbourhood.
}
\usage{
tsk(query, obs, subset = T, nmin = 3, nmax = 20, xcoord = "x", 
    ycoord = "y", tcoord = "t", zcoord = "z", vth = NULL, vlen = NULL, 
    llim = c(3, 3), verbose = T, Large = 2000, future = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{query}{
a data frame containing query spatiotemporal locations
}
  \item{obs}{
a data frame containing spatiotemporal locations and observed data
}
  \item{subset}{
logical; for local kriging; if TRUE only observations within the distances of estimated spatial and temporal sills from the prediction location are used for prediction}
  \item{nmin}{
  for local kriging: if the number of neighbors after subset is less than nmin, a missing value will be generated
  }
  \item{nmax}{
  for local kriging: the number of nearest observations that should be used for a kriging prediction, by default all observations are used.
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
	\item{future}{
	include observed points in future relative to query points.
}

}
\details{
Function implements global time space kriging based on a product sum model and support kriging in a local neighborhood.

If the number of neighbors exceeds a user-specified upper limit (\code{Large}), neighbors are sub-sampled in a balanced way to reduce the neighborhood size. 

Four variogram models: Gaussian, exponential, spherical and Matern are automatically fit to the empirical space and time variogram in the first lag. The range parameter is estimated from the first distance lag where the empirical variogram exceeds 80\% of the maximum. Weighted least square is then used to estimate the nugget and partial sill parameters. Model with minimal residual sum of squares between the empirical and fitted variogram is chosen as the variogram model. 

Field names for geographic coordinates and time stamps must match between query and observed data frames.
}
\value{
a list of a matrix \code{krig}, containing the prediction and prediction standard error and a flag denoting the reason for un-successful prediction:
\tabular{ll}{
0 \tab valid prediction\cr
1 \tab not enough temporal neighbors \cr
2 \tab not enough spatial neighbors \cr
3 \tab not enough neighbors \cr
4 \tab variogram could not be fit \cr
}
a list of estimated time space variogram, and a list of fitted parameter values of the product sum variogram model.

}
\references{
Iaco, S. De & Myers, D. E. & Posa, D., 2001. "Space-time analysis using a general product-sum model," Statistics & Probability Letters, Elsevier, vol. 52(1), pages 21-28, March.

Kumar, N., et al. (2013). "Satellite-based PM concentrations and their application to COPD in Cleveland, OH." Journal of Exposure Science and Environmental Epidemiology 23(6): 637-646.
	
Liang, D. and N. Kumar (2013). "Time-space Kriging to address the spatiotemporal misalignment in the large datasets." Atmospheric Environment 72: 60-69.
}
\author{
Naresh Kumar (NKumar@med.miami.edu)
Dong Liang (dliang@umces.edu)
}

\seealso{
\code{krigeST} in \pkg{gstat}
}

\examples{
## load the data
data(ex)
data(epa_cl)
## apply log transformation
obs[,'pr_pm25'] = log(obs[,'pr_pm25'])
## run kriging
system.time(out <- tsk(ex2.query[1:2,],obs,zcoord='pr_pm25',Large=400))
out$krig
}


