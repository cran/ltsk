lk <- function(query,obs,th,xcoord='x',ycoord='y',zcoord='z',vlen=15)
{
  seed <- round(runif(1) * 1000000)
  l.query <- check_input(query,xcoord,ycoord,'tt__tt',zcoord)
  l.obs <- check_input(obs,xcoord,ycoord,'tt__tt',zcoord)
  set.seed(seed=seed)
  th <- th^2 ## sqaured the threshold to accomodate ANN search
  
  ## setup the file names 
  ## l.dir <- tempdir()
  ## l.config <- paste(l.dir,"/config.txt",sep="")
  l.config <- "./config.txt"
  #l.query_file <- paste(l.dir, "/query.csv",sep="")
  #l.data_file <- paste(l.dir,"/data.csv",sep="")
  ff <- file(l.config,"w")
  ## l.wd <- getwd()
  ## setwd(l.dir)
  #write.table(x=l.query[,c(xcoord,ycoord,zcoord)],file=l.query_file,sep=",",row.names=F,col.names=F)
  #write.table(x=l.obs[,c(xcoord,ycoord,zcoord)],file=l.data_file,sep=",",row.names=F,col.names=F)
  #write the config.txt file under the current working directory
  cat("method=2D \n","query_file=./query.csv \n",
	"data_file=./data.csv \n",sep="",file=ff)
  cat("query_longitude_column=0\n","query_latitude_column=1\n",sep="",file=ff)
  cat("neighbour_longitude_column=0 \n","neighbour_latitude_column=1\n",sep="",file=ff)
  cat("more_than_one_z=0 \n",file=ff)
  cat("neighbour_z_colmn =2 \n",file=ff)
  cat("the_number_of_round =1 \n",file=ff)
  cat("dist_range=",th,"\n",file=ff)
  cat("dist_interval=",th,"\n",file=ff)
  cat("Rs=",vlen,"\n",file=ff)
  cat("Rt=",vlen,"\n",file=ff)
  cat("MAXPTS_ANN=",nrow(obs)+1000,"\n",file=ff)
  close(ff)
  ## write the data under the current working directory
  ## invoke the local kriging library
  .C("lk_main",
	as.double(l.obs[,xcoord]),
	as.integer(nrow(l.obs)),
	as.double(l.obs[,ycoord]),
	as.integer(nrow(l.obs)),
	as.double(l.obs[,zcoord]),
	as.integer(nrow(l.obs)),
	as.double(l.query[,xcoord]),
	as.integer(nrow(l.query)),
	as.double(l.query[,ycoord]),
	as.integer(nrow(l.query)))
   ## collect the result
  l.out_file <- paste('output',0,sep='_')
  l.out <- read.table(l.out_file,sep=",",na.strings="null")
  ## setwd(l.wd)
  colnames(l.out) <- c('null','krig','sigma','Hs','nugget','psill','model')
  l.out[,-1]
  ## garbage collection
}
