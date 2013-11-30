bootstrap_dyn <-function(function_,predFunction=NULL,getAdja_function, getBaseline_function,kfold=NULL,times,obs,n,b,K,delta,lambda,annot,annot_node,T_,active_mu,active_sd,inactive_mu,inactive_sd,prior=NULL,sourceNode=NULL,sinkNode=NULL,allint=FALSE,allpos=FALSE,muPgene=FALSE,muPgk=FALSE,muPgt=FALSE,muPgkt=FALSE,deltaPk=FALSE,deltaPt=FALSE,deltaPkt=FALSE)
{
	sq_err_all <- vector()

	edges_all = vector()
	baseline_all = vector()
	
  for(iter in 1:times)
  {
		#bootstrap replicates
		obs_modified = apply(obs, c(1,2,3), boot_mean)
		
		# solve LP problem
		res <- function_(obs=obs_modified,delta=delta,lambda=lambda,b=b,n=n,K=K,T_=T_,annot,prior=prior,sourceNode=sourceNode,sinkNode=sinkNode,all.int=allint,all.pos=allpos,deltaPk=deltaPk,deltaPt=deltaPt,deltaPkt=deltaPkt)

		#get adjacency matrix and baseline values
		adja <- getAdja_function(res$res,n,annot)
		baseline <- getBaseline_function(res$res,n=n)

		edges_all <- rbind(edges_all,c(t(adja)))
		baseline_all <- rbind(baseline_all, baseline)
		
		act_mat = calcActivation(adja, b, n, K)
		obs_inferred = getObsMat(act_mat,active_mu,active_sd,inactive_mu,inactive_sd)
		
		# calculate the squared error between the original obs matrix at t=T
		#  averaged over the replicates and the obs matrix generated from
		#  the inferred network at t=T 
		sq_err_all <- rbind(sq_err_all,c((obs_inferred-apply(obs[,,T_,], c(1,2), mean))^2))

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
