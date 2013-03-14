loocv <-function(function_,times,obs,n,b,K,delta,lambda,lambda2,annot,annot_node,T_=NULL,previousNet=NULL,baseline_=NULL,previousBaseline=NULL,active_mu,active_sd,inactive_mu,inactive_sd,prior=NULL,sourceNode=NULL,sinkNode=NULL,allint=FALSE,allpos=FALSE)
{
#	print("active_mu")
#	print(active_mu)
#	print("active_sd")
#	print(active_sd)
#	print("inactive_mu")
#	print(inactive_mu)
#	print("inactive_sd")
#	print(inactive_sd)
	
  kds <- matrix(b,nrow=n,ncol=K)
  # elements to leave out (each element at least once)
  looc <- cbind(rep(seq(1,dim(obs)[1]),dim(obs)[2]),rep(seq(1,dim(obs)[2]),rep(dim(obs)[1],dim(obs)[2])))
#  print("looc")
#  print(looc)
  adja_sum <- adja_num <- matrix(0,ncol=n,nrow=n)
  edges_all <- sq_err <- baseline_all <- vector()
  
  # observation of genes n after knockdowns k 
  for(x in 1:dim(looc)[1])
  {
		sq_err_tmp <- vector()
		obs_modified <- obs
		## randomly select an entry to be missing
		rem_gene <- looc[x,1] 	# which gene is removed
		rem_kd <- looc[x,2]		# in which knockdown
		obs_modified[rem_gene,rem_kd] <- NA
		ele <- obs[rem_gene,rem_kd]
		# mache nur, wenn der datenpunkt nicht eh schon NA ist
		if(!is.na(ele)){
			## do ILP
			res <- function_(obs=obs_modified,delta=delta,lambda=lambda,lambda2=lambda2,b=b,n=n,K=K,T_=T_,annot=annot,previousNet=previousNet,baseline=baseline_,previousBaseline=previousBaseline,prior=prior,sourceNode=sourceNode,sinkNode=sinkNode,all.int=allint,all.pos=allpos)
						
			adja <- getAdja(res,n)
			baseline <- getBaseline(res,n=n)
			for(i in 1:times){
				## calculate mean squared error of predicted and observed
				predict <- calcPredictionLOOCV(kds,adja,obs_modified,delta,rem_kd,rem_gene,active_mu,active_sd,inactive_mu,inactive_sd)
				sq_err_tmp <- c(sq_err_tmp,((predict-ele)^2))
			}
			res <- NA
			## calculate statistics on learned edges
			edges_all <- rbind(edges_all,c(t(adja)))
			baseline_all <- rbind(baseline_all, baseline)
#			print("loocv")
#			print(baseline)
#			print(baseline_all)
			sq_err <- c(sq_err,mean(sq_err_tmp,na.rm=T))
		}
  }
#   adja_mu <- adja_sum/(times*dim(looc)[1])
#   adja_prob <- adja_num/(times*dim(looc)[1])
  tmp1 <- rep(annot_node,rep(n,n))
  tmp2 <- rep(annot_node,n)
  id_selfloop <- which(tmp1==tmp2)
  tmp <- paste(tmp1,tmp2,sep="->")
  edges_all <- edges_all[,-id_selfloop]
  colnames(edges_all) <- tmp[-id_selfloop]
  MSE <- mean(sq_err,na.rm=T)
  
  return(list(MSE=MSE,edges_all=edges_all, baseline_all=baseline_all))
}







#-----------------------------------------------------------------------
# dyn version from 11/3/2013
#-----------------------------------------------------------------------


loocv_dyn <-function(function_,lambda2_,times,obs,n,b,K,delta,lambda,annot,annot_node,T_=NULL,previousNet=NULL,baseline_=NULL,previousBaseline=NULL,active_mu,active_sd,inactive_mu,inactive_sd,prior=NULL,sourceNode=NULL,sinkNode=NULL,allint=FALSE,allpos=FALSE)
{

	kds = array(NA, c(n,K,T_))
	for (t in 1:T_){
		kds[,,t] = matrix(b,nrow=n,ncol=K)
	}

  # elements to leave out (each element at least once)
  looc <- cbind(rep( seq(1,dim(obs)[1]), dim(obs)[2]),
								rep( seq(1,dim(obs)[2]), rep(dim(obs)[1],dim(obs)[2])),
								rep( seq(1,dim(obs)[3]), rep(dim(obs)[1]*dim(obs)[2],dim(obs)[3])))
								
#  print(looc)

  adja_sum <- adja_num <- matrix(0,ncol=n,nrow=n)
  edges_all <- sq_err <- baseline_all <- vector()
  
	# observation of genes n after knockdowns k 
	for(x in 1:dim(looc)[1])
	{
		sq_err_tmp <- vector()
		obs_modified <- obs
		## randomly select an entry to be missing
		rem_gene <- looc[x,1] 	# which gene is removed
		rem_kd <- looc[x,2]		# in which knockdown
		rem_t <- looc[x,3]	# in which time point
		obs_modified[rem_gene,rem_kd,rem_t] <- NA
		ele <- obs[rem_gene,rem_kd,rem_t]
		
		# mache nur, wenn der datenpunkt nicht eh schon NA ist
		if(!is.na(ele)){
			## do ILP
			res <- function_(obs_modified,delta,lambda=lambda,lambda2_,b,n,K,T_,annot,previousNet,baseline_,previousBaseline,prior=prior,sourceNode=sourceNode,sinkNode=sinkNode,all.int=allint,all.pos=allpos)
			
			adja <- getAdja(res,n)
			baseline <- getBaseline(res,n=n)
			
			for(i in 1:times){
				## calculate mean squared error of predicted and observed
				predict <- calcPredictionLOOCV_dyn(rem_t,kds,adja,obs_modified,delta,rem_kd,rem_gene,active_mu,active_sd,inactive_mu,inactive_sd)
				sq_err_tmp <- c(sq_err_tmp,((predict-ele)^2))
			}
			res <- NA
			## calculate statistics on learned edges
			edges_all <- rbind(edges_all,c(t(adja)))
			baseline_all <- rbind(baseline_all, baseline)
#			print("loocv")
#			print(baseline)
#			print(baseline_all)
			sq_err <- c(sq_err,mean(sq_err_tmp,na.rm=T))
		}
  }
#   adja_mu <- adja_sum/(times*dim(looc)[1])
#   adja_prob <- adja_num/(times*dim(looc)[1])
  tmp1 <- rep(annot_node,rep(n,n))
  tmp2 <- rep(annot_node,n)
  id_selfloop <- which(tmp1==tmp2)
  tmp <- paste(tmp1,tmp2,sep="->")
  edges_all <- edges_all[,-id_selfloop]
  colnames(edges_all) <- tmp[-id_selfloop]
  MSE <- mean(sq_err,na.rm=T)
  
  return(list(MSE=MSE,edges_all=edges_all, baseline_all=baseline_all))
}












#-----------------------------------------------------------------------
# nonIterative version
#-----------------------------------------------------------------------


loocv_nonIterative <-function(function_,lambda2_,times,obs,n,b,K,delta,lambda,annot,annot_node,T_=NULL,previousNet=NULL,baseline_=NULL,previousBaseline=NULL,active_mu,active_sd,inactive_mu,inactive_sd,prior=NULL,sourceNode=NULL,sinkNode=NULL,allint=FALSE,allpos=FALSE)
{

	kds = array(NA, c(n,K,T_))
	for (t in 1:T_){
		kds[,,t] = matrix(b,nrow=n,ncol=K)
	}

  # elements to leave out (each element at least once)
  looc <- cbind(rep( seq(1,dim(obs)[1]), dim(obs)[2]),
								rep( seq(1,dim(obs)[2]), rep(dim(obs)[1],dim(obs)[2])),
								rep( seq(1,dim(obs)[3]), rep(dim(obs)[1]*dim(obs)[2],dim(obs)[3])))
								
#  print(looc)

  adja_sum <- adja_num <- matrix(0,ncol=n,nrow=n)
  edges_all <- sq_err <- baseline_all <- vector()
  
  edges_all_list = list()
  baseline_all_list = list()
  
  for (t in 1:T_){
		edges_all_list[[t]] = matrix(c(NA), nrow=dim(looc)[1],ncol=n*n)
		baseline_all_list[[t]] = matrix(c(NA), nrow=dim(looc)[1],ncol=n)
  }
  
  
	# observation of genes n after knockdowns k 
	for(x in 1:dim(looc)[1])
	{
		sq_err_tmp <- vector()
		obs_modified <- obs
		## randomly select an entry to be missing
		rem_gene <- looc[x,1] 	# which gene is removed
		rem_kd <- looc[x,2]		# in which knockdown
		rem_t <- looc[x,3]	# in which time point
		obs_modified[rem_gene,rem_kd,rem_t] <- NA
		ele <- obs[rem_gene,rem_kd,rem_t]
#		print("rem")
#		print(rem_gene)
#		print(rem_kd)
#		print(rem_t)
		
		# mache nur, wenn der datenpunkt nicht eh schon NA ist
		if(!is.na(ele)){
			## do ILP
			res <- function_(obs_modified,delta,lambda=lambda,lambda2_,b,n,K,T_,annot,previousNet,baseline_,previousBaseline,prior=prior,sourceNode=sourceNode,sinkNode=sinkNode,all.int=allint,all.pos=allpos)

#			print(paste("lambda: ", lambda, sep=""))
			adja <- getAdja_nonIterative(res,n,T_)
			baseline <- getBaseline(res,n=n,T_=T_)
#			print(paste("x ", x, sep=""))
#			print(res$solution)
#			print("adja")
#			print(adja)
#			print("baseline")
#			print(baseline)
			
			for(i in 1:times){
				## calculate mean squared error of predicted and observed
				predict <- calcPredictionLOOCV_nonIterative(rem_t,kds,adja,obs_modified,delta,rem_kd,rem_gene,active_mu,active_sd,inactive_mu,inactive_sd)
				
				sq_err_tmp <- c(sq_err_tmp,((predict-ele)^2))
			}
			
#			print(paste("predict ", x, sep=""))
#			print(predict)
			res <- NA
			## calculate statistics on learned edges
			
			for ( t in 1:T_){
#				print(paste("transpose", t, sep=""))
#				print(c(t(adja[,,t])))
				edges_all_list[[t]][x,] <- c(t(adja[,,t]))
				baseline_all_list[[t]][x,] <- baseline
			}
			
#			print("edges_all_list")
#			print(edges_all_list)
#			print("baseline_all_list")
#			print(baseline_all_list)
#			print("loocv")
#			print(baseline)
#			print(baseline_all)
			sq_err <- c(sq_err,mean(sq_err_tmp,na.rm=T))
		}
	} # end x

#   adja_mu <- adja_sum/(times*dim(looc)[1])
#   adja_prob <- adja_num/(times*dim(looc)[1])
	
	
#	print("begin edges all")
#  print(edges_all_list)
  
  tmp1 <- rep(annot_node,rep(n,n))
  tmp2 <- rep(annot_node,n)
  id_selfloop <- which(tmp1==tmp2)
  tmp <- paste(tmp1,tmp2,sep="->")
  
  for (t in 1:T_){
		edges_all_list[[t]] <- edges_all_list[[t]][,-id_selfloop]
		colnames(edges_all_list[[t]]) <- tmp[-id_selfloop]
	}
	
#	print("final edges all")
#  print(edges_all_list)
  
  MSE <- mean(sq_err,na.rm=T)
  
  return(list(MSE=MSE,edges_all=edges_all_list, baseline_all=baseline_all_list))
}
