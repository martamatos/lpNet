runLP_dyn_Kfold_sevJobs_new = function(LPfunction, CVfunction, predFunction, kfold, geneState_, sd_all, totalruns, loocv_times, replnum, b_, n, K, T_, T_undNet, annot, annot_node, active_mu, inactive_mu, active_sd, inactive_sd, stepsize,outputDir,prior=NULL,startNode=NULL,endNode=NULL,allint=FALSE,allpos=FALSE,muPgene=FALSE,muPgk=FALSE,muPgt=FALSE,muPgkt=FALSE,deltaPk=FALSE,deltaPt=FALSE,deltaPkt=FALSE)
{

#--------------------------------------------------------
# for more than one time point
#--------------------------------------------------------
	print("T und net")
	print(T_undNet)
	act_mat = calcActivation(T_undNet,b_,n,K)
	
	sd_i = 1
	print("act_mat")
	print(act_mat)
	
	aucROC = list()
	aucPR = list()
	
	nw_all = list()
	baselines_all = list()
	
	obs = list()
	delta = list()

	
	
	for(std in sd_all)
	{
		baselines_all[[sd_i]] <- nw_all[[sd_i]] <- aucROC[[sd_i]] <- aucPR[[sd_i]] <- vector()
		
		edges = list()
		baselines = list()
		
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
					obs_tmp <- getObsMat(act_mat,active_mu,std,inactive_mu,std,geneState)
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
			ret <- doIt_dyn(LPfunction=LPfunction,CVfunction=CVfunction, predFunction=predFunction,loocv_times=loocv_times, kfold=kfold,stepsize=stepsize,obs=obs[[sd_i]],n=n,b=b,K=K,delta=delta[[sd_i]],annot=annot,annot_node=annot_node,T_=T_,prior=prior,startNode=startNode,endNode=endNode,allint=allint,active_mu=active_mu,active_sd=active_sd,inactive_mu=inactive_mu,inactive_sd=inactive_sd,muPgene=muPgene,muPgk=muPgk,muPgt=muPgt,muPgkt=muPgkt,deltaPk=deltaPk,deltaPt=deltaPt,deltaPkt=deltaPkt)

			edges[[i]] = ret$edges_all
			baselines[[i]] = ret$baseline_all
			
			colnames(T_undNet) <- rownames(T_undNet) <- annot_node

			## calculate Sensitivity and Specificity and ROC curve
			cv_roc <- calc_ROC(ret$edges_all,T_undNet,path="loocv_roc",triple=TRUE,plot=FALSE,t=t,std=std,outputDir=outputDir)
			aucROC[[sd_i]] <- c(aucROC[[sd_i]],cv_roc[[2]])
			aucPR[[sd_i]] <- c(aucPR[[sd_i]],cv_roc[[3]])

			
			# if returned best network is not null store it together with the baseline values
			if (!is.na(cv_roc$best_nw)){
				nw_all[[sd_i]] = rbind(nw_all[[sd_i]], cv_roc$best_nw)
				
				ret$baseline_all[which(ret$baseline_all==0)] = NA
				baselines_all[[sd_i]] = rbind(baselines_all[[sd_i]], apply(ret$baseline_all,2,mean,na.rm=T))
				baselines_all[[sd_i]][which(is.na(baselines_all[[sd_i]]))] = 0
			}

			
		} # end of total_runs
		
		save(edges,file=sprintf("%s/edges_std%s.Rdata", outputDir,std))
		save(baselines,file=sprintf("%s/baselines_std%s.Rdata", outputDir,std))
		
		sd_i = sd_i + 1
	} # end of sd_all
	
	return(list(aucROC=aucROC, aucPR=aucPR, nw_all=nw_all, baselines_all=baselines_all, obs=obs, delta=delta))
}
