runLP_dyn_Kfold_sevJobs_new = function(LPfunction, CVfunction, predFunction,getAdja_function, getBaseline_function, kfold, geneState_, sd_all, totalruns, loocv_times, replnum, b_, n, K, T_, T_undNet, annot, annot_node, active_mu, inactive_mu, active_sd, inactive_sd, stepsize,outputDir,prior=NULL,startNode=NULL,endNode=NULL,allint=FALSE,allpos=FALSE,mu_type=mu_type,delta_type=delta_type)
{

#--------------------------------------------------------
# for more than one time point
#--------------------------------------------------------
	print("T und net")
	print(T_undNet)
	act_mat = calcActivation_dyn(T_undNet,b_,n,K)
	
	sd_i = 1
	print("act_mat")
	print(act_mat)
	
	aucROC = list()
	aucPR = list()

	obs = list()
	delta = list()
	
	
	for(std in sd_all)
	{
		aucROC[[sd_i]] <- aucPR[[sd_i]] <- vector()
		
		edges = list()
		baselines = list()
		
		for(i in 1:totalruns)
		{
			# generate observation matrix and delta vector
			obs[[sd_i]] = array(NA, c(n,K,T_))
			tmp = c()
			
			for (t in 1:T_){
				geneState = geneState_[,,t]
				obs_all <-  vector()

				for(repl in 1:replnum)
				{
					obs_tmp <- getObsMat(act_mat,active_mu,std,inactive_mu,std,geneState)
					obs_all <- rbind(obs_all,c(obs_tmp))
					tmp <- cbind(tmp,rnorm(n,mean(c(active_mu,inactive_mu)),std))
				}
				obs[[sd_i]][,,t] <- matrix(apply(obs_all,2,mean,na.rm=T),nrow=n,ncol=K)
			}
			delta[[sd_i]] <- apply(tmp,1,mean,na.rm=T)

			
			print("obs")
			print(obs[[sd_i]] )
			print(delta[[sd_i]] )
			print(b_)
			print(geneState_)

			# run nonIterative model
			ret <- doIt_dyn(LPfunction=LPfunction,CVfunction=CVfunction, predFunction=predFunction,getAdja_function=getAdja_function, getBaseline_function=getBaseline_function,loocv_times=loocv_times, kfold=kfold,stepsize=stepsize,obs=obs[[sd_i]],n=n,b=b_,K=K,delta=delta[[sd_i]],annot=annot,annot_node=annot_node,T_=T_,prior=prior,startNode=startNode,endNode=endNode,allint=allint,active_mu=active_mu,active_sd=active_sd,inactive_mu=inactive_mu,inactive_sd=inactive_sd,mu_type=mu_type,delta_type=delta_type)

			edges[[i]] = ret$edges_all
			baselines[[i]] = ret$baseline_all
			
#			colnames(T_undNet) <- rownames(T_undNet) <- annot_node

#			# calculate Sensitivity and Specificity and ROC curve
#			cv_roc <- calc_ROC(edges_all=ret$edges_all,T_nw=T_undNet,path=NULL,triple=TRUE,plot=FALSE,sampleMAD=TRUE, method1=median, method2=mad, method2Times=1, septype="->")
			
#			aucROC[[sd_i]] <- c(aucROC[[sd_i]],cv_roc[[2]])
#			aucPR[[sd_i]] <- c(aucPR[[sd_i]],cv_roc[[3]])

			
		} # end of total_runs
		
		save(edges,file=sprintf("%s/edges_std%s.Rdata", outputDir,std))
		save(baselines,file=sprintf("%s/baselines_std%s.Rdata", outputDir,std))
		
		sd_i = sd_i + 1
	} # end of sd_all
	
	return(list(aucROC=aucROC, aucPR=aucPR, obs=obs, delta=delta))
}
