\name{lk}
\alias{lk}
\title{
local space kriging
}
\description{

This function implements ANN neighbor search, automatic variogram estimation
and ordinary Kriging to predict spatial data at multiple query locations.

}
\usage{
lk(query, obs, th, xcoord = "x", ycoord = "y", zcoord = "z",  vlen = 15)
}
\arguments{
  \item{query}{
data frame containing query spatial locations
}
  \item{obs}{
data frame containing spatial locations and observed data
}
  \item{th}{
a priori chosen distance thresholds for neighbor search
}
  \item{xcoord}{
field name for x coordinate in both query and obs
}
  \item{ycoord}{
field name for y coordinate in both query and obs}

  \item{zcoord}{
field name for data in obs}

  \item{vlen}{
 number of bins to compute local spatial variogram}
}
\details{

With the development of GIS and remote sensing, data sets are available over large spatial domain.
Ordinary kriging which involves the O(n^3) matrix operation, becomes computationally in-feasible for prediction. In addition to the computational issues, the stationarity assumption underlying ordinary kriging is usually violated for environmental data over large domain. 

In addition to Bayesian spatial modeling, local method offers a computational efficient approach for spatial prediction. This function implements a local Kriging method for spatial prediction from large data sets. At each location where prediction is required (termed as query point), a plausible and concervative distance threshold is used to select neighboring locations with observed spatial data. A local variogram is estimated to enable ordinary kriging within this neighborhood.

Implementation of the local Kriging requires specifying a distance threshold (th). It is an upper bound within which spatial auto-correlation is believed to exist. The distance threshold may be obtained from variogram estimates plus knowledge of the underlying spatial process. For the given threhold, ANN tree is used to efficiently identify neighbors.

Remote sensed data are dense within a small distance. Dense data lead to computational issue for variogram estimation and kriging even for local methods. In senarios when the number of neighbors exceeds a user-specified upper limit (Large), neighbors are sub-sampled in a spatially balanced way to reduce the neighborhood size. Sub-sampling leads to improved computational efficiency without loss of predictive accuracy due to the  strong spatial autocorrelation within short distances. 

In remote sensed applications, cloud cover could lead to missing data. This poses challenges in the subsequent variogram estimation and Kriging. Without further knowledge of the spatial process, query points within the area with sparse data are excluded. This function does not perform local kriging when the neighborhood contains less than \emph{five} distinct spatial locations with observed data in the neighborhood. But we recommend at least 20 data points for varioram estimation and Kriging.

Local empirical variogram is calculated using method of moments within the neighborhood. Four variogram models: Gaussian, exponential, spherical and Matern are fit to the empirical variogram. The range parameter is estimated from the first distance lag where the empirical variogram exceeds 80% of the maximum. Weighted least square is then used to estimate the nugget and partial sill parameters. Model with minimal residaul sum of sqaures between the empirical and fitted variogram is chosen as the variogram model. 

Ordinary kriging is used to predict the spatial process at each query point. The estimated range parameter is used to reduce the neighborhood. Only initial neighbors that are spatially correlated with the process at the query point are used in prediction. In case the estimated range is small, indicating a locally weak spatial process, a sample of 10% of the neighbor was used instead.
}
\value{
Kriging mean and standard error and quality flags.\tabular{ll}{
krig \tab predictive mean\cr
sigma \tab predictive standard error \cr
Hs \tab estimated range \cr
nugget \tab estimated nugget \cr
psill \tab estimated partial sill \cr
model \tab variogram model \cr
}
where the variogram model are
\tabular{ll}{
1 \tab Spherical \cr
2 \tab Gaussian \cr
3 \tab Exponential \cr
4 \tab Matern \cr
}
}
\examples{
data(epa_cl)
data(ex)
out <- lk(ex1.grid,ex1.data,0.3)
quilt.plot(ex1.grid$x,ex1.grid$y,out$krig,nrow=20,ncol=20,main='predicted PM[2.5] 01/01/2005')
}
\author{
Jun chen (wdidwlia@gmail.com)
Jin Chen (jc.chenjin@gmail.com)
Dong Liang (dong-liang@uiowa.edu)
}


