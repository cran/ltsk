dnb <-
function(query, obs,th,future=TRUE)
{
  ## query : coordinate and time stamp of the query point
  ## obs   : coordinates and time stamps of observed points
  ## th    : space and time thresholds
  ## value : coordinates and time stamps of potential neighbors
  if (is.data.frame(query)){
	  query <- as.matrix(query)
	  query <- as.vector(query)
  }
  if (is.data.frame(obs)){
  	obs <- as.matrix(obs)
  }
  if(future){
    tvec <- abs( query[3] - obs[,3] )
    iit <- which( tvec < th[2] )    
  }
  else{  
    tvec <- obs[,3]- query[3] 
    iit <- which( tvec > - th[2] & tvec <= 0 )
  }
  loc0 <- matrix(query[1:2], ncol=2)
  locs <- obs[iit,1:2]
  dvec <- as.vector(rdist(loc0,locs))
  iid <- which(dvec < th[1])
  iit[iid]
}
