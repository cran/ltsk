check_input <-
function(in_data,xcoord,ycoord,tcoord,zcoord)
{
   l.coords <- c(xcoord,ycoord,tcoord,zcoord)
   if( is.null(dim(in_data)))
   {
	    ## data is a vector
	    l.data <- matrix(in_data,ncol=length(in_data),nrow=1)
	    colnames(l.data) <- names(in_data)
   }
   else
   {
	    ## query is a matrix
	    l.data <- in_data
   }
   
   for (l.coord in l.coords)
   {
	  if (is.na(match(l.coord,colnames(l.data)))){
		  l.txt <- colnames(l.data)
		  l.data <- cbind(l.data,-99999)
		  ##warning(l.coord,"not found in input data!\n")
		  colnames(l.data) <- c(l.txt,l.coord)
	  }
   }		
   r <- as.matrix(l.data[,l.coords])
   
   ## check for duplicates
   isdup <- duplicated(r)
   if(mean(isdup)>0.1){
     arg <- deparse(substitute(in_data))
     warning(" over 10% duplicates in ",arg,", please check. \n")
   }
#    if(any(isdup)){
#     #if(mean(isdup)>0.5){ ## over 50% duplicated
# 	  arg <- deparse(substitute(in_data))
# 	  warning(sum(isdup), " duplicates removed from ", arg,"\n")
# 	  r <- r[!isdup,]
#    }
#    
   ## add an attribute that is removed if true
   attr(r,"remove.action") <- isdup
   r
}

check_na <- function(in_data,type_txt)
{
  if( is.null(dim(in_data)))
  {
    ## data is a vector
    l.data <- matrix(in_data,ncol=length(in_data),nrow=1)
    colnames(l.data) <- names(in_data)
  }
  else{
    l.data <- in_data
  }
  
  r <- na.omit(l.data)
  if(is.null(attr(l.data,"remove.action"))){
    ## no duplicate removed
    attr(r,"remove.action") <- rep(F,nrow(r))
  }else{
    ## some duplicate removed
    if(is.null(attr(r,"na.action"))){
      ## no missing value removed
      ## no action
    }else{
      ## missing value removed
      ii <- which(!attr(l.data,"remove.action")) ## current observations
      jj <- ii[attr(r,"na.action")]
      attr(l.data,"remove.action")[jj] <- T ## remove missing values
    }
    attr(r,"remove.action") <- attr(l.data,"remove.action")
  }
  r
  ## no_na <- na.omit(l.data)
  
  #dd <- attr(no_na,"na.action")
  #if (length(dd)>0){
	#	cat(length(dd)," ",type_txt," location removed.\n")
	#}
	## no_na
}