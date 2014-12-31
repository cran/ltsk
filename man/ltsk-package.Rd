\name{ltsk-package}
\alias{ltsk-package}
\docType{package}
\title{
Local Time Space Kriging
}
\description{
ltsk library is a collection of programs for implementing local spatial and local spatiotemporal Kriging. 
Unlike global Kriging, ltsk subsets the sample around a given location and time where predicted is needed; 
estimates variogram using the subset of sample data. 
Product-sum model is implemented and automatically estimated using the data points within the local neighbourhood. 
A unique advantage of ltsk is that it addresses non-stationarity,
 which is difficult to handle in large spatiotemporal dataset. 
}
\details{
\tabular{ll}{
Package: \tab ltsk\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-12-31\cr
License: \tab GPL2.0\cr
}
}
\author{
Naresh Kumar (NKumar@med.miami.edu)
Dong Liang (dliang@umces.edu)
Jun chen (wdidwlia@gmail.com)
Jin Chen (jc.chenjin@gmail.com)
}
\references{
Haas, Timothy C. "Local prediction of a spatio-temporal process with an application to wet sulfate deposition." Journal of the American Statistical Association 90.432 (1995): 1189-1199.

Iaco, S. De & Myers, D. E. & Posa, D., 2001. "Space-time analysis using a general product-sum model," Statistics & Probability Letters, Elsevier, vol. 52(1), pages 21-28, March.

Kumar, N., et al. (2013). "Satellite-based PM concentrations and their application to COPD in Cleveland, OH." Journal of Exposure Science and Environmental Epidemiology 23(6): 637-646.
	
Liang, D. and N. Kumar (2013). "Time-space Kriging to address the spatiotemporal misalignment in the large datasets." Atmospheric Environment 72: 60-69.
	
}
\keyword{package}

