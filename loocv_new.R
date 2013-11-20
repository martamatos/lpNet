loocv_LP <-
function(function_,predFunction,kfold=NULL,times,obs,n,b,K,delta,lambda,annot,annot_node,T_=NULL,active_mu,active_sd,inactive_mu,inactive_sd,prior=NULL,sourceNode=NULL,sinkNode=NULL,allint=FALSE,allpos=FALSE,muPgene=FALSE,muPgk=FALSE,muPgt=FALSE,muPgkt=FALSE,deltaPk=FALSE,deltaPt=FALSE,deltaPkt=FALSE)
{
	
  # elements to leave out (each element at least once)
  looc <- cbind(rep(seq(1,dim(obs)[1]),dim(obs)[2]),rep(seq(1,dim(obs)[2]),rep(dim(obs)[1],dim(obs)[2])))

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
			res <- function_(obs=obs_modified,delta=delta,lambda=lambda,b=b,n=n,K=K,T_=T_,annot,prior=prior,sourceNode=sourceNode,sinkNode=sinkNode,all.int=allint,all.pos=allpos,deltaPk=deltaPk,deltaPt=deltaPt,deltaPkt=deltaPkt)
						
			adja <- getAdja(res,n)
			baseline <- getBaseline(res,n=n)
			res <- NA
			
			for(i in 1:times){
				## calculate mean squared error of predicted and observed
				predict <- predFunction(b=b,n=n,K=K,adja=adja,baseline=baseline,obs=obs_modified,delta=delta,rem_gene=rem_gene, rem_k=rem_kd,active_mu=active_mu,active_sd=active_sd,inactive_mu=inactive_mu,inactive_sd=inactive_sd,muPgene=muPgene,muPgk=muPgk)
				sq_err_tmp <- c(sq_err_tmp,((predict-ele)^2))
			}
			
			## calculate statistics on learned edges
			edges_all <- rbind(edges_all,c(t(adja)))
			baseline_all <- rbind(baseline_all, baseline)

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




loocv_dyn <-function(function_,predFunction,getAdja_function, getBaseline_function,kfold=NULL,times,obs,n,b,K,delta,lambda,annot,annot_node,T_,active_mu,active_sd,inactive_mu,inactive_sd,prior=NULL,sourceNode=NULL,sinkNode=NULL,allint=FALSE,allpos=FALSE,muPgene=FALSE,muPgk=FALSE,muPgt=FALSE,muPgkt=FALSE,deltaPk=FALSE,deltaPt=FALSE,deltaPkt=FALSE)
{

  # elements to leave out (each element at least once)
  looc <- cbind(rep( seq(1,dim(obs)[1]), dim(obs)[2]),
								rep( seq(1,dim(obs)[2]), rep(dim(obs)[1],dim(obs)[2])),
								rep( seq(1,dim(obs)[3]), rep(dim(obs)[1]*dim(obs)[2],dim(obs)[3])))
								
  adja_sum <- adja_num <- matrix(0,ncol=n,nrow=n)
  edges_all <- sq_err <- baseline_all <- vector()
	
	# dont remove observations from first time point
  looc = looc[-which(looc[,3]==1, arr.ind=T),]

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
			res <- function_(obs=obs_modified,delta=delta,lambda=lambda,b=b,n=n,K=K,T_=T_,annot,prior=prior,sourceNode=sourceNode,sinkNode=sinkNode,all.int=allint,all.pos=allpos,deltaPk=deltaPk,deltaPt=deltaPt,deltaPkt=deltaPkt)

			adja <- getAdja_function(res$res,n)
			baseline <- getBaseline_function(res$res,n=n)
			res <- NA
			
			for(i in 1:times){
				## calculate mean squared error of predicted and observed
				predict <- predFunction(b=b,n=n,K=K,adja=adja,baseline=baseline,obs=obs_modified,delta=delta,rem_gene=rem_gene, rem_k=rem_kd, rem_t=rem_t,active_mu=active_mu,active_sd=active_sd,inactive_mu=inactive_mu,inactive_sd=inactive_sd,muPgene=muPgene,muPgk=muPgk,muPgt=muPgt,muPgkt=muPgkt)
				sq_err_tmp <- c(sq_err_tmp,((predict-ele)^2))
			}
			
			## calculate statistics on learned edges
			edges_all <- rbind(edges_all,c(t(adja)))
			baseline_all <- rbind(baseline_all, baseline)
			sq_err <- c(sq_err,mean(sq_err_tmp,na.rm=T))
		}
	}

  tmp1 <- rep(annot_node,rep(n,n))
  tmp2 <- rep(annot_node,n)
  id_selfloop <- which(tmp1==tmp2)
  tmp <- paste(tmp1,tmp2,sep="->")
  edges_all <- edges_all[,-id_selfloop]
  colnames(edges_all) <- tmp[-id_selfloop]
  MSE <- mean(sq_err,na.rm=T)

  return(list(MSE=MSE,edges_all=edges_all, baseline_all=baseline_all))
}
