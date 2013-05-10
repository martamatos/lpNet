runLP_dyn_Kfold_sevJobs_2lev = function(function_, lambda2_, kfold, geneState_, sd_all, totalruns, loocv_times, replnum, b_, n, K, T_, T_nw_,T_undNet, annot, annot_node, active_mu, inactive_mu, inactivePknock_mu, active_sd, inactive_sd, inactivePknock_sd, graphSize_,outputDir)
{

#--------------------------------------------------------
# for more than one time point
#--------------------------------------------------------
	print("T und net")
	print(T_undNet)
	act_mat = calcActivation(T_undNet,b_,n,K)
	
	sd_i = 1
	print(act_mat)
	
	aucROC = list()
	aucROC_noNA = list()
	aucPR = list()
	
	randROC_list = list()
	randPR_list = list()
	
	nw_all = list()
	baselines_all = list()
	
	obs = list()
	delta = list()

	
	
	for(std in sd_all)
	{
		baselines_all[[sd_i]] <- nw_all[[sd_i]] <- aucROC[[sd_i]] <- aucPR[[sd_i]] <- aucROC_noNA[[sd_i]] <- nw_size <- vector()
		
		randROC_list[[sd_i]] <- vector()
		randPR_list[[sd_i]] <- vector()
		edges = list()
		
		for(i in 1:totalruns)
		{
			# generate observation matrix and delta vector
			obs[[sd_i]] = array(NA, c(n,K,T_))
			tmp = c()
			
			for (t in 1:T_){
				b = b_
				geneState = geneState_[,,t]
				obs_all <-  vector()

				for(repl in 1:replnum)
				{
					obs_tmp <- getObsMat_2lev(act_mat,active_mu,std,inactive_mu,std,inactivePknock_mu,std,geneState)
					obs_all <- rbind(obs_all,c(obs_tmp))
					tmp <- cbind(tmp,rnorm(n,mean(c(active_mu,inactive_mu)),std))
				}
				obs[[sd_i]][,,t] <- matrix(apply(obs_all,2,mean,na.rm=T),nrow=n,ncol=K)
			}
			delta[[sd_i]] <- apply(tmp,1,mean,na.rm=T)

			
			print("obs")
			print(obs)
			print(delta)
			print(b)
			print(geneState_)

			# run nonIterative model
			ret <- doIt_dyn_kfold_2lev(function_,kfold,loocv_times,obs[[sd_i]],n,b,K,delta[[sd_i]],lambda=lamd,lambda2_,annot,annot_node,T_,previousNet,baseline,previousBaseline,prior=NULL,startNode=NULL,endNode=NULL,allint=FALSE,active_mu,active_sd,inactive_mu,inactive_sd, inactivePknock_mu=inactivePknock_mu, inactivePknock_sd=inactivePknock_sd)

			edges[[i]] = ret$edges_all
			
#				print("ret$edges_all")
#				print(ret$edges_all)
			
			colnames(T_undNet) <- rownames(T_undNet) <- annot_node


			## calculate Sensitivity and Specificity and ROC curve
			loocv_roc <- calc_ROC(ret$edges_all,T_undNet,path="loocv_roc",triple=TRUE,plot=FALSE,t=t,std=std,outputDir=outputDir)

			aucROC[[sd_i]] <- c(aucROC[[sd_i]],loocv_roc[[2]])
			aucPR[[sd_i]] <- c(aucPR[[sd_i]],loocv_roc[[3]])

			
			# if returned best network is not null store it together with the baseline values
			if (!is.na(loocv_roc$best_nw)){
				nw_all[[sd_i]] = rbind(nw_all[[sd_i]], loocv_roc$best_nw)
				
				ret$baseline_all[which(ret$baseline_all==0)] = NA
				baselines_all[[sd_i]] = rbind(baselines_all[[sd_i]], apply(ret$baseline_all,2,mean,na.rm=T))
				baselines_all[[sd_i]][which(is.nan(baselines_all[[sd_i]]))] = 0

				aucROC_noNA[[sd_i]] = c(aucROC_noNA[[sd_i]],loocv_roc[[2]])
			}

			
		} # end of total_runs
		
		save(edges,file=sprintf("%s/edges_std%s.Rdata", outputDir,std))
		
		sd_i = sd_i + 1
	} # end of sd_all
	
	return(list(randROC_list=randROC_list, randPR_list=randPR_list, aucROC=aucROC, aucPR=aucPR, aucROC_noNA=aucROC_noNA, nw_all=nw_all, baselines_all=baselines_all, obs=obs, delta=delta))
}
