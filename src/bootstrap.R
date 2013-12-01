bootstrap_dyn <-function(function_,predFunction=NULL,getAdja_function, getBaseline_function,kfold=NULL,times,obs,n,b,K,delta,lambda,annot,annot_node,T_,active_mu,active_sd,inactive_mu,inactive_sd,prior=NULL,sourceNode=NULL,sinkNode=NULL,allint=FALSE,allpos=FALSE,mu_type,delta_type)
{
	
	sq_err_all <- vector()

	edges_all = vector()
	baseline_all = vector()
	
  for(iter in 1:times)
  {
		#bootstrap replicates
		obs_modified = apply(obs, c(1,2,3), boot_mean)
		
		# solve LP problem
		res <- function_(obs=obs_modified,delta=delta,lambda=lambda,b=b,n=n,K=K,T_=T_,annot,prior=prior,sourceNode=sourceNode,sinkNode=sinkNode,all.int=allint,all.pos=allpos,delta_type=delta_type)

		#get adjacency matrix and baseline values
		adja <- getAdja_function(res$res,n,annot)
		baseline <- getBaseline_function(res$res,n=n)
		edges_all <- rbind(edges_all,c(t(adja)))
		baseline_all <- rbind(baseline_all, baseline)
		
		act_mat = calcActivation(adja, b, n, K)
		obs_inferred = getObsMat(act_mat,active_mu,active_sd,inactive_mu,inactive_sd, mu_type)
		
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

bootstrap_dyn_ddepn<-function(function_,predFunction=NULL,getAdja_function, getBaseline_function,kfold=NULL,times,obs,n,b,K,delta,lambda,annot,annot_node,T_,active_mu,active_sd,inactive_mu,inactive_sd,prior=NULL,sourceNode=NULL,sinkNode=NULL,allint=FALSE,allpos=FALSE,mu_type,delta_type)
{
	delta = NULL
	active_mu = NULL
	active_sd = NULL
	inactive_mu = NULL
	inactive_sd = NULL
	
	sq_err_all <- vector()

	edges_all = vector()
	baseline_all = vector()
	
  for(iter in 1:times)
  {
		#bootstrap replicates
		obs_modified = apply(obs, c(1,2,3), boot_mean)
		obs_modified = round(obs_modified, digits=5)
		
		delta = obs_modified[,,1]
		delta[3:14,]=delta[3:14,]+0.00001
		delta = round(delta, digits=5)

		active_mu = matrix(NA,nrow=16, ncol=3)
		active_sd = matrix(0.01,nrow=16, ncol=3)
		inactive_mu = matrix(NA,nrow=16, ncol=3)
		inactive_sd = matrix(0.01,nrow=16, ncol=3)

		for (e in 1:K){
			for (i in 1:n){ 
				vec = obs_modified[i,e,]
				active_mu[i,e] =  mean(vec[which(vec >= delta[i,e])], na.rm=T)
				inactive_mu[i,e] = mean(vec[which(vec < delta[i,e])], na.rm=T)
				
			}
		}
		
		for (i in 1:n){
			active_mu[i,which(is.na(active_mu[i,]))]=mean(active_mu[i,], na.rm=T)
			inactive_mu[i,which(is.na(inactive_mu[i,]))]=mean(inactive_mu[i,], na.rm=T)
		}
		
		# solve LP problem
		res <- function_(obs=obs_modified,delta=delta,lambda=lambda,b=b,n=n,K=K,T_=T_,annot,prior=prior,sourceNode=sourceNode,sinkNode=sinkNode,all.int=allint,all.pos=allpos,delta_type=delta_type)

		#get adjacency matrix and baseline values
		adja <- getAdja_function(res$res,n,annot)
		baseline <- getBaseline_function(res$res,n=n)
		edges_all <- rbind(edges_all,c(t(adja)))
		baseline_all <- rbind(baseline_all, baseline)
		
		act_mat = calcActivation(adja, b, n, K)
		obs_inferred = getObsMat(act_mat,active_mu,active_sd,inactive_mu,inactive_sd, mu_type)
		
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
