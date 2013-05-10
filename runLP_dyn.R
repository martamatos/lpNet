runLP_dyn = function(function_, lambda2_, geneState_, sd_all, totalruns, loocv_times, replnum, b_, n, K, T_, T_nw_,T_undNet, annot, annot_node, active_mu, inactive_mu, active_sd, inactive_sd,graphSize_,outputDir)
{

#--------------------------------------------------------
# for more than one time point
#--------------------------------------------------------
	print("T und net")
	print(T_undNet)
	act_mat = calcActivation(T_undNet,b_,n,K)
	sd_i = 1
	print(act_mat)
	auc_list = matrix(c(NA), nrow=length(sd_all), ncol=4)
	auc_list_median  = matrix(c(NA), nrow=length(sd_all), ncol=10)
	randROC_list = vector()
	randPR_list = vector()
	
	
	
	for (k in 1:K){
		system(sprintf('mkdir %s/k%s', outputDir,k))
	}

	for (std in sd_all){
		for (k in 1:K){
			system(sprintf('mkdir %s/k%s/std%s', outputDir,k,std))
		}
	}

	for(std in sd_all)
	{
		baselines_all <- nw_all <- aucROC <- aucPR <- aucROC_noNA <- nw_size <- vector()
		
		edges = list()
		
		for(i in 1:totalruns)
		{
			# generate observation matrix and delta vector
			obs = array(NA, c(n,K,T_))
			tmp = c()
			
			for (t in 1:T_){
				b = b_[t,]
				geneState = geneState_[,,t]
				obs_all <-  vector()

				for(repl in 1:replnum)
				{
					obs_tmp <- getObsMat(act_mat,active_mu,std,inactive_mu,std,geneState)
					obs_all <- rbind(obs_all,c(obs_tmp))
					tmp <- cbind(tmp,rnorm(n,mean(c(active_mu,inactive_mu)),std))
				}
				obs[,,t] <- matrix(apply(obs_all,2,mean,na.rm=T),nrow=n,ncol=K)
			}
			delta <- apply(tmp,1,mean,na.rm=T)
			b = b_
			
			print("obs")
			print(obs)
			print(delta)
			print(b)
			print(geneState_)

			# run nonIterative model
			ret <- doIt_dyn(function_,lambda2_,loocv_times,obs,n,b,K,delta,lambda=lamd,annot,annot_node,T_,previousNet,baseline,previousBaseline,prior=NULL,startNode=NULL,endNode=NULL,allint=FALSE,active_mu,active_sd,inactive_mu,inactive_sd)
			
				print("ret$edges_all")
				print(ret$edges_all)
			edges[[i]] = ret$edges_all
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
		
		save(edges,file=sprintf("%s/edges_std%s.Rdata", outputDir,std))
		
		print(aucROC)
		print(aucPR)
		
		# get the median/mean values for AUC-ROC AUC-PR
		auc_list[sd_i,] <- c(mean(aucROC,na.rm=T),mean(aucPR,na.rm=T),sd(aucROC,na.rm=T),sd(aucPR,na.rm=T))
		auc_list_median[sd_i,] <- c(quantile(aucROC,na.rm=T),quantile(aucPR,na.rm=T))

		
		
		for (t in 1:T_){
			for (k in 1:K){
				bestNW_dyn(n,obs[,k,t], delta, nw_all,baselines_all, std, aucROC, aucROC_noNA, aucPR, t, k,graphSize_,outputDir, function_)
			}
		}
		
		sd_i = sd_i + 1
	} # end of sd_all
	
	# print AUC results 
#	print(randROC_list)
	printAUC_und(auc_list,sd_all,randROC_list,randPR_list,totalruns,t,outputDir)
	printAUC_median_und(auc_list_median,sd_all,randROC_list,randPR_list,totalruns,t,outputDir)
	
}
