bootstrap_dyn <-function(function_,predFunction=NULL,getAdja_function, getBaseline_function,kfold=NULL,times,obs,n,b,K,delta,lambda,annot,annot_node,T_,active_mu=NULL,active_sd=NULL,inactive_mu=NULL,inactive_sd=NULL,prior=NULL,sourceNode=NULL,sinkNode=NULL,allint=FALSE,allpos=FALSE,muPgene=FALSE,muPgk=FALSE,muPgt=FALSE,muPgkt=FALSE,deltaPk=FALSE,deltaPt=FALSE,deltaPkt=FALSE)
{
	sq_err_all <- vector()

	edges_all = vector()
	baseline_all = vector()
	
  for(iter in 1:times)
  {
		#bootstrap replicates
		obs_modified = apply(obs, c(1,2,3), boot_mean)
		
		## solve LP problem
		res <- function_(obs=obs_modified,delta=delta,lambda=lambda,b=b,n=n,K=K,T_=T_,annot,prior=prior,sourceNode=sourceNode,sinkNode=sinkNode,all.int=allint,all.pos=allpos,deltaPk=deltaPk,deltaPt=deltaPt,deltaPkt=deltaPkt)

		#get adjacency matrix and baseline values
		adja <- getAdja_function(res$res,n,annot)
		baseline <- getBaseline_function(res$res,n=n)

		edges_all <- rbind(edges_all,c(t(adja)))
		baseline_all <- rbind(baseline_all, baseline)
	
		# calculate the squared error between the bootstraped observation matrix 
		#  and the original obs matrix averaged over the replicates
		sq_err_all <- rbind(sq_err_all,c((obs_modified-apply(obs, c(1,2,3), mean))^2))
  } # end iter

  sq_err <- apply(sq_err_all,2,mean,na.rm=T)
  tmp1 <- rep(annot_node,rep(n,n))
  tmp2 <- rep(annot_node,n)
  id_selfloop <- which(tmp1==tmp2)
  tmp <- paste(tmp1,tmp2,sep="->")
  edges_all <- edges_all[,-id_selfloop]
  colnames(edges_all) <- tmp[-id_selfloop]
  MSE <- mean(sq_err,na.rm=T)

  return(list(MSE=MSE,edges_all=edges_all,baseline_all=baseline_all))
}

# small function to bootstrap vec
boot_mean = function(vec){
    mean(sample(vec, size=length(vec), replace=TRUE))
}
