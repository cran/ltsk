\name{dnb}
\alias{dnb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Time space neighbor search
}
\description{
Identify locations adjacent in both space and time to a given location.
}
\usage{
dnb(query, obs, th)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{query}{
	a vector of x/y coordinates and time stamp of the query point
}
  \item{obs}{
	coordinates and time stamps of the observed points
}
  \item{th}{
	distance and time thresholds 
}
}
\details{
Let query point denote the X/Y coordinates and time stamp for which neighbour search is needed.  Query point may represent the spatio-temporal location of a case. Neighbouring data are needed for exposure assessment.

Let observed points denote all locations (X/Y coordinates and time stamps) with observed exposure data. 
  
This function identified the observed data within a given distance and time thresholds. Observed data within a given radius and time period both prior to and after the query location will be included.
}
\value{
A vector of locations of observed data that are within the given distance and time threshold.
}
\author{
Dong Liang (dong-liang@uiowa.edu)
}
\note{
This is a naive neighbour search. It works well for data that are dense in space but sparse in time. Future work involves range tree implementation for large and dense spatiotemporal data.
}


\examples{
data(epa_cl)
coords <- c('x','y','t')
ii <- dnb(query[1,coords],obs[,coords],c(0.1,10))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
