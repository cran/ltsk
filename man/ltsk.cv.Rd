\name{ltsk.cv}
\alias{ltsk.cv}
\title{
Local Time and Space Kriging Cross Validation, n-Fold or Leave-one-out
}
\description{
Cross validation functions for local time space kriging}
\usage{
ltsk.cv(nfold, obs, th, nbins, part=NULL,zcoord = "z",...)
}
\arguments{
  \item{nfold}{
integer, apply n-fold cross validation; if larger than number of observed data, apply leave-one-out cross validation}
  \item{obs}{
data frame containing spatiotemporal locations and observed data
}
  \item{th}{
vector of length two; a priori chosen distance threshold and time lag for neighbor search
}
  \item{nbins}{
vector of length two; a priori chosen bins to divide distance threshold and time lag equally}
  \item{part}{
vector of random digits between 1 and \code{nfold}; if NULL, it was sampled with replacement from \code{seq(1,nfold)} of length \code{nrow(obs)}
}
  \item{zcoord}{
	character constant, the field name for data in \code{obs}
}
  \item{...}{
other arguments that will be passed to \code{cltsk}
}

}
\details{
Leave-one-out cross validation visits a data point, and predicts the value at that location by leaving out the observed value, and proceeds with the next data point. N-fold cross validation makes a partitions the data set in N parts. For all observations in a part, predictions are made based on the remaining N-1 parts; this is repeated for each of the N parts.}
\value{
a matrix of the cross validation residual, each column corresponds to a given distance threshold and time lag; a data frame containing the summary statistics of the cross validation residuals, including number of non-missing kriging, the sum of square prediction errors and the mean square prediction errors. Each individual row is a combination of distance threshold and time lag.
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

\examples{
## load the data
set.seed(123)
data(epa_cl)
ii= with(obs,which(amonth==5 & aday <13)) ## first week of Januray 2005;
x=obs[sample(ii,400),]
## apply log transformation
x[,'pr_pm25'] = log(x[,'pr_pm25'])
## run kriging
out <- ltsk.cv(nfold=10,obs=x,th=c(0.10,10),nbins=c(2,2),zcoord='pr_pm25',verbose=FALSE)
}


