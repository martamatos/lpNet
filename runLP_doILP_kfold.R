runLP_doILP_kfold = function(function_, lambda2_, kfold, geneState_, sd_all, totalruns, loocv_times, replnum, b_, n, K, T_, T_nw_,T_undNet, annot, annot_node, active_mu, inactive_mu, active_sd, inactive_sd,graphSize_,outputDir)
{

#--------------------------------------------------------
# for more than one time point
#--------------------------------------------------------
	print("T und net")
	print(T_undNet)
	act_mat = calcActivation(T_undNet,b_,n,K)
	sd_i = 1

	auc_list = matrix(c(NA), nrow=length(sd_all), ncol=4)
	auc_list_median  = matrix(c(NA), nrow=length(sd_all), ncol=10)
	randROC_list = vector()
	randPR_list = vector()

	
	for(std in sd_all)
	{
		baselines_all <- nw_all <- aucROC <- aucPR <- aucROC_noNA <- nw_size <- vector()
		
		for(i in 1:totalruns)
		{
			# generate observation matrix and delta vector
			obs = matrix(NA, nrow=n, ncol=K)
			tmp = c()
			
			b = b_
			geneState = geneState_
			
			obs_all <-  vector()

			for(repl in 1:replnum)
			{
				obs_tmp <- getObsMat(act_mat,active_mu,std,inactive_mu,std,geneState)
				obs_all <- rbind(obs_all,c(obs_tmp))
				tmp <- cbind(tmp,rnorm(n,mean(c(active_mu,inactive_mu)),std))
			}
			obs <- matrix(apply(obs_all,2,mean,na.rm=T),nrow=n,ncol=K)
			
			delta <- apply(tmp,1,mean,na.rm=T)
			
			print("obs")
			print(obs)
			print(delta)
			print(b)
			print(geneState)

			# run nonIterative model
			ret <- doIt_kfold(function_=doILP,kfold=kfold,loocv_times=loocv_times,obs=obs,n=n,b=b,K=K,delta=delta,lambda=lamd,lambda2=lambda2,annot=annot,annot_node=annot_node,T_=T_,previousNet=previousNet,baseline=baseline,previousBaseline=previousBaseline,prior=NULL,startNode=NULL,endNode=NULL,allint=FALSE,active_mu=active_mu,active_sd=active_sd,inactive_mu=inactive_mu,inactive_sd=inactive_sd)
			 
			
#				print("ret$edges_all")
#				print(ret$edges_all)
			
			colnames(T_undNet) <- rownames(T_undNet) <- annot_node
			randRoc = calcRandROC(T_undNet,100)
			randROC_list <- c(randROC_list, randRoc[[1]])
			randPR_list <- c(randPR_list, randRoc[[2]])

			## calculate Sensitivity and Specificity and ROC curve
			loocv_roc <- calc_ROC(ret$edges_all,T_undNet,path="loocv_roc",triple=TRUE,plot=FALSE,t=t,std=std,outputDir=outputDir)
			nw_size <- c(nw_size,sum(loocv_roc[[4]]!=0))
			aucROC <- c(aucROC,loocv_roc[[2]])
			aucPR <- c(aucPR,loocv_roc[[3]])

			
			# if returned best network is not null store it together with the baseline values
			if (!is.na(loocv_roc$best_nw)){
				nw_all = rbind(nw_all, loocv_roc$best_nw)
				
				ret$baseline_all[which(ret$baseline_all==0)] = NA
				baselines_all = rbind(baselines_all, apply(ret$baseline_all,2,mean,na.rm=T))
				baselines_all[which(is.nan(baselines_all))] = 0

				aucROC_noNA = c(aucROC_noNA,loocv_roc[[2]])
			}

			
		} # end of total_runs
		
		
		# get the median/mean values for AUC-ROC AUC-PR
		auc_list[sd_i,] <- c(mean(aucROC,na.rm=T),mean(aucPR,na.rm=T),sd(aucROC,na.rm=T),sd(aucPR,na.rm=T))
		auc_list_median[sd_i,] <- c(quantile(aucROC,na.rm=T),quantile(aucPR,na.rm=T))
		
		

		system(sprintf('mkdir %s/std%s', outputDir,std))

		
		
		for (t in 1:T_){
			bestNW(n, nw_all, baselines_all, std, aucROC, aucROC_noNA, aucPR, t, graphSize_, outputDir, function_)
		}
		
		sd_i = sd_i + 1
	} # end of sd_all
	
	# print AUC results 
#	print(randROC_list)
	printAUC_und(auc_list,sd_all,randROC_list,randPR_list,totalruns,t,outputDir)
	printAUC_median_und(auc_list_median,sd_all,randROC_list,randPR_list,totalruns,t,outputDir)
	
}
