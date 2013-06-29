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
		l.data <- cbind(l.data,NA)
		warning(l.coord,"not found in input data!\n")
		colnames(l.data) <- c(l.txt,l.coord)
	}
   }		
   as.matrix(l.data[,l.coords])
}
