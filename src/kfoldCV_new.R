kfoldCV_LP <-
function(function_,predFunction,kfold,times,obs,n,b,K,delta,lambda,annot,annot_node,T_=NULL,active_mu,active_sd,inactive_mu,inactive_sd,prior=NULL,sourceNode=NULL,sinkNode=NULL,allint=FALSE,allpos=FALSE,mu_type,delta_type)
{
  # define k-fold groups: stratified
  obs_kfold <- list()
  num <- (dim(obs)[1]*dim(obs)[2])
  le <- ceiling(num/kfold)
  sq_err_all <- vector()
  edges_all <- baseline_all <- vector()
  sq_err <- vector()

	
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
			res <- function_(obs=obs_modified,delta=delta,lambda=lambda,b=b,n=n,K=K,T_=T_,annot,prior=prior,sourceNode=sourceNode,sinkNode=sinkNode,all.int=allint,all.pos=allpos,delta_type=delta_type)
		
			adja <- getAdja(res,n)
			baseline <- getBaseline(res,n=n)
			## calculate statistics on learned edges
			edges_all <- rbind(edges_all,c(t(adja)))
			baseline_all <- rbind(baseline_all, baseline)
			
			## calculate mean squared error of predicted and observed
			predict <- predFunction(b=b,n=n,K=K,adja=adja,baseline=baseline,obs=obs_modified,delta=delta,rem_entries=rem_entries,rem_entries_vec=rem_entries_vec,active_mu=active_mu,active_sd=active_sd,inactive_mu=inactive_mu,inactive_sd=inactive_sd,mu_type=mu_type)
			
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

  
  return(list(MSE=MSE,edges_all=edges_all, baseline_all=baseline_all))
}




kfoldCV_dyn <-function(function_,predFunction,getAdja_function, getBaseline_function,kfold,times,obs,n,b,K,delta,lambda,annot,annot_node,T_,active_mu,active_sd,inactive_mu,inactive_sd,prior=NULL,sourceNode=NULL,sinkNode=NULL,allint=FALSE,allpos=FALSE,mu_type,delta_type)
{

  # define k-fold groups: stratified
  obs_kfold <- list()
  num <- (dim(obs)[1]*dim(obs)[2]*(dim(obs)[3]-1)) # observations from t=1 are not removed
  le <- ceiling(num/kfold)
  sq_err_all <- vector()
  sq_err <- vector()

	edges_all = vector()
	baseline_all = vector()
	
	# trick so that observations from t=1 are not removed, part 1
	obs_complete = obs
	obst1 = obs[,,1]
	obs = obs[,,-1]
	
	NAentries = which(is.na(obs))

	
  for(k in 1:times){
		tmps <- sample(seq(1,kfold),kfold)
		sq_err_tmp <- vector()
		
		for(j in 1:kfold){
			tmp <- c(tmps[j:kfold],tmps[1:j-1])
			obs_order <- order(obs)
			obs_kfold[[j]] <- array(NA,c(n,K,T_-1))
			for(i in 1:le){
				if(num>=i){
					obs_kfold[[j]][obs_order[tmp[1]]] <- obs[obs_order[tmp[1]]]
				}
			# 	  print(tmp)
				obs_order <- obs_order[-tmp]
				tmp <- c(tmp[-1],tmp[1])
			}
			obs_kfold[[j]] <- array(obs_kfold[[j]],c(n,K,T_-1))
		}
		
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
			
			train_data <- array(train_data,c(n,K,T_-1))
			obs_modified <- train_data
		
			rem_entries = which(is.na(obs_modified), arr.ind=TRUE)
			rem_entries_vec = which(is.na(obs_modified))

			rowsToRemove = which(rem_entries_vec %in% intersect(rem_entries_vec, NAentries))

			
			if (length(rowsToRemove)>0){
				rem_entries = rem_entries[-rowsToRemove,]
				rem_entries_vec = rem_entries_vec[-rowsToRemove]
			}
			rem_entries[,3] = rem_entries[,3] + 1
			rem_entries_vec = rem_entries_vec + n*K

			
			# trick so that observations from t=1 are not removed, part 2
			obstemp = array(NA, c(n,K,T_))
			obstemp[,,1] = obst1
			obstemp[,,2:T_] = obs_modified
			obs_modified = obstemp

			
			## do ILP
			res <- function_(obs=obs_modified,delta=delta,lambda=lambda,b=b,n=n,K=K,T_=T_,annot,prior=prior,sourceNode=sourceNode,sinkNode=sinkNode,all.int=allint,all.pos=allpos,delta_type=delta_type)

			adja <- getAdja_function(res$res,n)
			baseline <- getBaseline_function(res$res,n=n)
			## calculate statistics on learned edges
			edges_all <- rbind(edges_all,c(t(adja)))
			baseline_all <- rbind(baseline_all, baseline)
			
			## calculate mean squared error of predicted and observed
			predict <- predFunction(b=b,n=n,K=K,adja=adja,baseline=baseline,obs=obs_modified,delta=delta,rem_entries=rem_entries,rem_entries_vec=rem_entries_vec,active_mu=active_mu,active_sd=active_sd,inactive_mu=inactive_mu,inactive_sd=inactive_sd,mu_type=mu_type)

#			ids_rem <- which(is.na(obs_modified))
			ids_rem <- rem_entries_vec
			sq_err_tmp <- c(sq_err_tmp,((predict[ids_rem]-obs_complete[ids_rem])^2))
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
  
  return(list(MSE=MSE,edges_all=edges_all,baseline_all=baseline_all))
}

