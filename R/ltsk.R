ltsk <-
function(query,obs,th,xcoord='x',ycoord='y',tcoord='t',zcoord='z',
		vth=NULL,vlen=NULL,llim=c(3,3),verbose=T,Large=2000,nproc=NULL)
{
 seed <- round(runif(1) * 1000000)
 l.query <- check_input(query,xcoord,ycoord,tcoord,zcoord)
 l.obs <- check_input(obs,xcoord,ycoord,tcoord,zcoord)
 if(is.null(nproc)){
	set.seed(seed=seed)
	## Use single core mode
 	out <- apply(l.query,1,working.ltsk,obs=l.obs,th=th,vth=vth,vlen=vlen,
		llim=llim,verbose=verbose,Large=Large)
 	out <- t(out)
 }
 else{
	## Use multiple core mode
 	cl <- makeSOCKcluster(nproc)
	clusterSetupRNG(cl)
	clusterSetupRNGstream(cl,seed=rep(seed,6))
 	pwd <- getwd()
 	clusterCall(cl,setwd,dir=pwd)
	clusterEvalQ(cl,library(ltsk))
 	#nmlist <- list('l.obs','th','xcoord','ycoord','tcoord','zcoord',
	#		'vth','vlen','llim','verbose','Large','check_input')
	#for(var in nmlist)
	#	assign(var,get(var), envir = .GlobalEnv)		
 	#clusterExport(cl,nmlist) ## Export only get variables in GlobalEnv
	#clusterEvalQ(cl,check_input(obs=obs,xcoord=xcoord,ycoord=ycoord,
	#		tcoord=tcoord,zcoord=zcoord))
 	#out <- parRapply(cl=cl,x=l.query,fun=working.ltsk,
	#		obs=l.obs,th=th,vth=vth,vlen=vlen,
	#		llim=llim,verbose=verbose,Large=Large)
	ll.query <- as.list(data.frame(t(l.query)))
	ll.obs <- vector('list',length(ll.query))
	for(i in 1:length(ll.query)){
		ii <- dnb(ll.query[[i]],l.obs,th)
		ll.obs[[i]] <- l.obs[ii,]
	}
	ll.args <- list(vth=vth,vlen=vlen,llim=llim,verbose=verbose,Large=Large)
	out1 <- clusterMap(cl=cl,fun=working.ltsk.par,ll.query,ll.obs,MoreArgs=ll.args)
 	stopCluster(cl)
 	#out <- matrix(out,ncol=3,byrow=T)
	out <- matrix(unlist(out1),ncol=3,byrow=T)
 }
 colnames(out) <- c('fit','se','flag')
 cbind(query,out)
}
