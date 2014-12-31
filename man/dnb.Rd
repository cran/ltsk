\name{dnb}
\alias{dnb}
\title{
Search Neighbours in Time and Space Within Specified Ranges
}
\description{
A brute force neighbor search implementation to identify observed data points within a given distance around location and time interval.
}
\usage{
dnb(query, obs, th, future=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{query}{
	a vector; the x, y coordinates and the time stamp of the query point
}
  \item{obs}{
	a matrix; the x, y coordinates and time stamps of the spatiotemporal locations
}
  \item{th}{
	a vector; the distance threshold and time lag
}
	\item{future}{
	logical, whether include observed spatiotemporal points future in time relative to the query spatiotemporal location.
}
}
\details{

Implementation involves first calculating the time lags between query point and observed data (with locational coordinates and time); for observed locations within time lag of query, the function calculates the Euclidean distances between query location and all potential neighbors and select those within specified distance threshold.

The future argument can be used to exclude data in the future in neighbor search. This is useful in an extrapolation application.
}
\value{
A vector, row numbers in the observed data matrix, that are within the given distance threshold and time lag of the query location.
}
\author{
Dong Liang (dliang@umces.edu)
}
\note{
For large dataset, use ANN (for spatial kriging) and Range Tree for spatiotemporal Kriging.
}

\seealso{
\code{get.knn} in \pkg{FNN}
}

\examples{
data(epa_cl)
coords <- c('x','y','t')
ii <- dnb(query[1,coords],obs[,coords],c(0.1,10))
}
