firstpeak <-
function(dist,gamma)
{
 ## dist: distance intervals
 ## gamma : estimated semi-variogram ( with missing values)
 ## value: first peak of gamma
 ii <- which(is.na(gamma))
 if( length(ii) > 0){
	f <- approxfun(dist, gamma,rule=2)
	gamma[ii] <- f(dist[ii]) 
 }
 maxg <- max(gamma)
 min(which( gamma > (maxg *.8)))
}
