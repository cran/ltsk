\name{lk}
\alias{lk}
\title{
Local Ordinary Kriging
}
\description{

Function implements ANN based neighbor search, automatic variogram estimation and ordinary Kriging at query locations.

}
\usage{
lk(query, obs, th, xcoord = "x", ycoord = "y", zcoord = "z",  vlen = 15, cl=NULL)
}
\arguments{
  \item{query}{
a data frame containing query spatial locations
}
  \item{obs}{
a data frame containing spatial locations and observed data
}
  \item{th}{
a vector, distance threshold for neighbor search
}
  \item{xcoord}{
a character constant, field name for x coordinate in both \code{query} and \code{obs}
}
  \item{ycoord}{
a character constant, field name for y coordinate in both \code{query} and \code{obs}}

  \item{zcoord}{
a character constant, field name for data in \code{obs} }

  \item{vlen}{
 number of bins to compute local spatial variogram}

 \item{cl}{
  a cluster object for multi-core computing (default number of cores in the PC minue one), 0 means single core
  }
 }
\details{

This function implements a local spatial Kriging method for spatial prediction using large data sets. At each location where prediction is required (termed as query point), a conservative distance threshold is used to select neighboring locations from the observed (or sample) data, which has XY coordinates.

A local variogram is used to estimate ordinary Kriging within this neighborhood.

Implementation of the local Kriging requires specifying a distance threshold (th). It is an upper bound within which spatial auto-correlation is believed to exist. The distance threshold may be obtained from variogram estimates or chosen based on prior knowledge of the underlying spatial process. For the given threshold, ANN tree is used to efficiently identify neighbors.

This function does not perform local kriging when the neighborhood contains less than \emph{five} distinct spatial locations with observed data in the neighborhood. More data points are recommended for robust varioram estimation.

Four variogram models: Gaussian, exponential, spherical and Matern are automatically fit to the empirical variogram. The range parameter is estimated from the first distance lag where the empirical variogram exceeds 80% of the maximum. Weighted least square is then used to estimate the nugget and partial sill parameters. Model with minimal residual sum of squares between the empirical and fitted variogram is chosen as the variogram model. 

The estimated range parameter is used to reduce the neighborhood. Only initial neighbors that are spatially correlated with the process at the query point are used in prediction. In case the estimated range is small, indicating a locally weak spatial process, a sample of 10% of the neighbor was used instead.

Field names for geographic coordinates must match between query and observed data frames.
}
\value{
kriging prediction and standard deviation and quality flags.\tabular{ll}{
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

\seealso{
  \code{autoKrige} in \pkg{automap}
}

\examples{
data(epa_cl)
data(ex)
out <- lk(ex1.grid[1:10,],ex1.data,0.5,cl=0)
}
\author{
Jun chen (wdidwlia@gmail.com)
Jin Chen (jc.chenjin@gmail.com)
Naresh Kumar (NKumar@med.miami.edu)
Dong Liang (dliang@umces.edu)
}


