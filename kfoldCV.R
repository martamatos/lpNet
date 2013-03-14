kfoldCV <-
function(function_,lambda2_,kfold,times,obs,n,b,K,delta,lambda,annot,annot_node,T_=NULL,previousNet=NULL,baseline_=NULL,previousBaseline=NULL,active_mu,active_sd,inactive_mu,inactive_sd,prior=NULL,sourceNode=NULL,sinkNode=NULL,allint=FALSE,allpos=FALSE)
{
  kds <- matrix(b,nrow=n,ncol=K)
  # define k-fold groups: stratified
  obs_kfold <- list()
  num <- (dim(obs)[1]*dim(obs)[2])
  le <- ceiling(num/kfold)
  sq_err_all <- vector()
  edges_all <- baseline_all <- vector()
  sq_err <- vector()
  
#  print("kfold")
#  print(num)
#  print(kfold)
#  print(le)
  
  for(k in 1:times){
		tmps <- sample(seq(1,kfold),kfold)
		sq_err_tmp <- vector()
		
		for(j in 1:kfold){
			tmp <- c(tmps[j:kfold],tmps[1:j-1])
			obs_order <- order(obs)
			obs_kfold[[j]] <- matrix(NA,nrow=n,ncol=K)
			
			for(i in 1:le){
				if(num>=i){
					obs_kfold[[j]][obs_order[tmp[1]]] <- obs[obs_order[tmp[1]]]
				}
			# 	  print(tmp)
				obs_order <- obs_order[-tmp]
				tmp <- c(tmp[-1],tmp[1])
			}
			obs_kfold[[j]] <- matrix(obs_kfold[[j]],nrow=n,ncol=K)
		}
		
		#   mean(obs_kfold[[1]],na.rm=T)
		## make crossvalidation
		adja_sum <- adja_num <- matrix(0,ncol=n,nrow=n)
		
		
		for(x in 1:kfold){
			test_ids <- seq(1,kfold)[-x]
			train_tmp <- vector()
			
			for(i in test_ids){
				train_tmp <- rbind(train_tmp,c(obs_kfold[[i]]))
			}
			train_data <- rep(NA,dim(train_tmp)[2])
			
			for(i in 1:dim(train_tmp)[2]){
				if(!all(is.na(train_tmp[,i]))){
					train_data[i] <- sum(train_tmp[,i],na.rm=T)
				}
			}
			train_data <- matrix(train_data,nrow=n,ncol=K)
			obs_modified <- train_data
		
			## do ILP
			res <- function_(obs_modified,delta,lambda=lambda,lambda2_,b,n,K,T_,annot,previousNet,baseline_,previousBaseline,prior=prior,sourceNode=sourceNode,sinkNode=sinkNode,all.int=allint,all.pos=allpos)
		
			adja <- getAdja(res,n)
			baseline <- getBaseline(res,n=n,T_=NULL)
			## calculate statistics on learned edges
			edges_all <- rbind(edges_all,c(t(adja)))
			baseline_all <- rbind(baseline_all, baseline)
			
			## calculate mean squared error of predicted and observed
			predict <- calcPredictionKfoldCV(adja,b,n,K,active_mu,active_sd,inactive_mu,inactive_sd)
			ids_rem <- which(is.na(obs_modified))
			sq_err_tmp <- c(sq_err_tmp,((predict[ids_rem]-obs[ids_rem])^2))
		}
		sq_err_all <- rbind(sq_err_all,sq_err_tmp)
  } # end times
  sq_err <- apply(sq_err_all,2,mean,na.rm=T)
  tmp1 <- rep(annot_node,rep(n,n))
  tmp2 <- rep(annot_node,n)
  id_selfloop <- which(tmp1==tmp2)
  tmp <- paste(tmp1,tmp2,sep="->")
  edges_all <- edges_all[,-id_selfloop]
  colnames(edges_all) <- tmp[-id_selfloop]
  MSE <- mean(sq_err,na.rm=T)
  
#  print("baseline_all")
#  print(baseline_all)
#  print(edges_all)
  
  return(list(MSE=MSE,edges_all=edges_all, baseline_all=baseline_all))
}







kfoldCV_dyn <-
function(function_,lambda2_,kfold,times,obs,n,b,K,delta,lambda,annot,annot_node,T_=NULL,previousNet=NULL,baseline_=NULL,previousBaseline=NULL,active_mu,active_sd,inactive_mu,inactive_sd,prior=NULL,sourceNode=NULL,sinkNode=NULL,allint=FALSE,allpos=FALSE)
{

	kds <- matrix(b,nrow=n,ncol=K)
  # define k-fold groups: stratified
  obs_kfold <- list()
  num <- (dim(obs)[1]*dim(obs)[2]*dim(obs)[3])
  le <- ceiling(num/kfold)
  sq_err_all <- vector()
  sq_err <- vector()
#  print(obs)

	edges_all = vector()
	baseline_all = vector()
	
	
#	print("cenas")
#	print(le)
#	print(obs)
  
  for(k in 1:times){
		tmps <- sample(seq(1,kfold),kfold)
		sq_err_tmp <- vector()
		
		for(j in 1:kfold){
			tmp <- c(tmps[j:kfold],tmps[1:j-1])
			obs_order <- order(obs)
			obs_kfold[[j]] <- array(NA,c(n,K,T_))
			for(i in 1:le){
				if(num>=i){
					obs_kfold[[j]][obs_order[tmp[1]]] <- obs[obs_order[tmp[1]]]
				}
			# 	  print(tmp)
				obs_order <- obs_order[-tmp]
				tmp <- c(tmp[-1],tmp[1])
			}
			obs_kfold[[j]] <- array(obs_kfold[[j]],c(n,K,T_))
		}
		
		#   mean(obs_kfold[[1]],na.rm=T)
		## make crossvalidation
		adja_sum <- adja_num <- matrix(0,ncol=n,nrow=n)
		
		
		for(x in 1:kfold){
			test_ids <- seq(1,kfold)[-x]
			train_tmp <- vector()
			
			for(i in test_ids){
				train_tmp <- rbind(train_tmp,c(obs_kfold[[i]]))
			}
			train_data <- rep(NA,dim(train_tmp)[2])
			
			for(i in 1:dim(train_tmp)[2]){
				if(!all(is.na(train_tmp[,i]))){
					train_data[i] <- sum(train_tmp[,i],na.rm=T)
				}
			}
			
			train_data <- array(train_data,c(n,K,T_))
			obs_modified <- train_data
#			print(paste("kfold: ", x, sep=""))
#			print(paste("time ", k, sep=""))
#			print("obs_modified")
#			print(obs_modified)
#			print(obs)
			
			## do ILP
			res <- function_(obs_modified,delta,lambda=lambda,lambda2_,b,n,K,T_,annot,previousNet,baseline_,previousBaseline,prior=prior,sourceNode=sourceNode,sinkNode=sinkNode,all.int=allint,all.pos=allpos)
			
			adja <- getAdja(res,n)
#			print("adja")
#			print(adja)
			baseline <- getBaseline(res,n=n,T_=NULL)
			## calculate statistics on learned edges
			edges_all <- rbind(edges_all,c(t(adja)))
			baseline_all <- rbind(baseline_all, baseline)
			
			## calculate mean squared error of predicted and observed
			predict <- calcPredictionKfoldCV_dyn(T_,adja,b,n,K,active_mu,active_sd,inactive_mu,inactive_sd)
			ids_rem <- which(is.na(obs_modified))
			sq_err_tmp <- c(sq_err_tmp,((predict[ids_rem]-obs[ids_rem])^2))
		}
		sq_err_all <- rbind(sq_err_all,sq_err_tmp)
  } # end times
  sq_err <- apply(sq_err_all,2,mean,na.rm=T)
  tmp1 <- rep(annot_node,rep(n,n))
  tmp2 <- rep(annot_node,n)
  id_selfloop <- which(tmp1==tmp2)
  tmp <- paste(tmp1,tmp2,sep="->")
  edges_all <- edges_all[,-id_selfloop]
  colnames(edges_all) <- tmp[-id_selfloop]
  MSE <- mean(sq_err,na.rm=T)
  
#  print("edges_all")
#  print(edges_all)
  
  return(list(MSE=MSE,edges_all=edges_all,baseline_all=baseline_all))
}











kfoldCV_nonIterative <-
#(times,obs,n,b,K,delta,lambda,annot,annot_node,kfold,active_mu,active_sd,
#         inactive_mu,inactive_sd,prior=NULL,sourceNode=NULL,sinkNode=NULL,allint=FALSE,allpos=FALSE)
function(function_,lambda2_,kfold,times,obs,n,b,K,delta,lambda,annot,annot_node,T_=NULL,previousNet=NULL,baseline_=NULL,previousBaseline=NULL,active_mu,active_sd,inactive_mu,inactive_sd,prior=NULL,sourceNode=NULL,sinkNode=NULL,allint=FALSE,allpos=FALSE)
{

	kds <- matrix(b,nrow=n,ncol=K)
  # define k-fold groups: stratified
  obs_kfold <- list()
  num <- (dim(obs)[1]*dim(obs)[2]*dim(obs)[3])
  le <- ceiling(num/kfold)
  sq_err_all <- vector()
  sq_err <- vector()
#  print(obs)

	edges_all_list = list()
	baseline_all_list = list()
	for (t in 1:T_){
		edges_all_list[[t]] = matrix(c(NA), nrow=kfold,ncol=n*n)
		baseline_all_list[[t]] = matrix(c(NA), nrow=kfold,ncol=n)
	}
	
#	print("cenas")
#	print(le)
#	print(obs)
  
  for(k in 1:times){
		tmps <- sample(seq(1,kfold),kfold)
		sq_err_tmp <- vector()
		
		for(j in 1:kfold){
			tmp <- c(tmps[j:kfold],tmps[1:j-1])
			obs_order <- order(obs)
			obs_kfold[[j]] <- array(NA,c(n,K,T_))
			for(i in 1:le){
				if(num>=i){
					obs_kfold[[j]][obs_order[tmp[1]]] <- obs[obs_order[tmp[1]]]
				}
			# 	  print(tmp)
				obs_order <- obs_order[-tmp]
				tmp <- c(tmp[-1],tmp[1])
			}
			obs_kfold[[j]] <- array(obs_kfold[[j]],c(n,K,T_))
		}
		
		#   mean(obs_kfold[[1]],na.rm=T)
		## make crossvalidation
		adja_sum <- adja_num <- matrix(0,ncol=n,nrow=n)
		
		
		for(x in 1:kfold){
			test_ids <- seq(1,kfold)[-x]
			train_tmp <- vector()
			
			for(i in test_ids){
				train_tmp <- rbind(train_tmp,c(obs_kfold[[i]]))
			}
			train_data <- rep(NA,dim(train_tmp)[2])
			
			for(i in 1:dim(train_tmp)[2]){
				if(!all(is.na(train_tmp[,i]))){
					train_data[i] <- sum(train_tmp[,i],na.rm=T)
				}
			}
			
			train_data <- array(train_data,c(n,K,T_))
			obs_modified <- train_data
#			print(paste("kfold: ", x, sep=""))
#			print(paste("time ", k, sep=""))
#			print("obs_modified")
#			print(obs_modified)
#			print(obs)
			
			## do ILP
			res <- function_(obs_modified,delta,lambda=lambda,lambda2_,b,n,K,T_,annot,previousNet,baseline_,previousBaseline,prior=prior,sourceNode=sourceNode,sinkNode=sinkNode,all.int=allint,all.pos=allpos)
			
			adja <- getAdja_nonIterative(res,n,T_)
			baseline <- getBaseline(res,n=n,T_=T_)
			
#			print("adja")
#			print(adja)
			for ( t in 1:T_){
#				print(paste("transpose", t, sep=""))
#				print(c(t(adja[,,t])))
				edges_all_list[[t]][x,] <- c(t(adja[,,t]))
				baseline_all_list[[t]][x,] <- baseline
			}
			
			
			## calculate statistics on learned edges
#			edges_all <- rbind(edges_all,c(t(adja)))
			
			## calculate mean squared error of predicted and observed
			predict <- calcPredictionKfoldCV_nonIterative(T_,adja,b,n,K,active_mu,active_sd,inactive_mu,inactive_sd)
			
			ids_rem <- which(is.na(obs_modified))
#			print("ids_rem")
#			print(ids_rem)
#			print(obs[ids_rem])
			sq_err_tmp <- c(sq_err_tmp,((predict[ids_rem]-obs[ids_rem])^2))
		} # end kfold
		sq_err_all <- rbind(sq_err_all,sq_err_tmp)
		
  } # end times
  
  sq_err <- apply(sq_err_all,2,mean,na.rm=T)
  tmp1 <- rep(annot_node,rep(n,n))
  tmp2 <- rep(annot_node,n)
  id_selfloop <- which(tmp1==tmp2)
  tmp <- paste(tmp1,tmp2,sep="->")
  
#  print("edges_all_list")
#  print(edges_all_list)
  for (t in 1:T_){
		edges_all_list[[t]] <- edges_all_list[[t]][,-id_selfloop]
		colnames(edges_all_list[[t]]) <- tmp[-id_selfloop]
	}
#	print(edges_all_list)
  
  MSE <- mean(sq_err,na.rm=T)
#  stop("cenas")
  return(list(MSE=MSE,edges_all=edges_all_list,baseline_all=baseline_all_list))
}
