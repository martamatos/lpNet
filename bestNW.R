options(digits=2) # so that results are printed only with 2 decimal places

bestNW = function(n, nw_all, baselines_all, std, aucROC, aucROC_noNA, aucPR, i, graphSize_, outputDir, function_)
{
	# print best network for each value of std
	print("-------------------------------------------------------------")
	print(sprintf('t: %s, sd: %s outputDir: %s',i,std, outputDir))
	print(aucROC)
#	print(aucROC_noNA)
	
	# safeguard for when the program finds only empty nets
	if (length(aucROC_noNA) == 0){
		nw_final = matrix(0, nrow=n, ncol=n)
		baseline_final = matrix(0, nrow=1, ncol=n)
		
		write.matrix(nw_final, file = sprintf("%s/adj_t=%s_std=%s",outputDir,i, std))
		write.matrix(baseline_final, file = sprintf("%s/baseline_t=%s_std=%s",outputDir,i, std))
		
		
		system(sprintf('python net5.py  printGraph_%s %s/adj_t=%s_std=%s %s/std%s/graph%s_%s',graphSize_, outputDir,i,std, outputDir,std,i,std), intern=TRUE)
		
		print(paste("max auc (zero): ", 0 , sep=""))
		print(nw_final)
		print(baseline_final)
	}
	else{
	
		ind_baseline = which(aucROC_noNA==max(aucROC_noNA,na.rm = TRUE))
		ind01 = n * (which(aucROC_noNA==max(aucROC_noNA,na.rm = TRUE))-1) + 1
		ind02 = ind01 + n -1
		nw_mean = baseline_mean= c()
		
		
		try(
			if (which(aucROC==max(aucROC,na.rm = TRUE)) != which(aucPR==max(aucPR,na.rm = TRUE))){
				print(" nw for best aucROC != nw for best aucPR")
			}, FALSE)
		
		
		for ( j in 1:length(ind01) ){
			nw_tmp = nw_all[ind01[j]:ind02[j],]
			nw_mean = rbind(nw_mean,c(t(nw_tmp)))
			baseline_tmp = baselines_all[ind_baseline[j],]
			baseline_mean = rbind(baseline_mean, baseline_tmp)
		}
		
	#	print(nw_all)
	#	print(nw_mean)
		nw_final = matrix(apply(nw_mean,2,median,na.rm=T),nrow=n,ncol=n,byrow=T)
		baseline_final = apply(baseline_mean, 2, median, na.rm=T)

		write.matrix(nw_final, file = sprintf("%s/adj_t=%s_std=%s",outputDir,i, std))
		write.matrix(baseline_final, file = sprintf("%s/baseline_t=%s_std=%s",outputDir,i, std))
		
		
		system(sprintf('python net5.py  printGraph_%s %s/adj_t=%s_std=%s %s/std%s/graph%s_%s', graphSize_, outputDir,i,std, outputDir,std,i,std), intern=TRUE)
		
		print(paste("max auc: ",max(aucROC_noNA,na.rm = TRUE), sep=""))
		print(nw_final)
		print(baseline_final)
	}
	return(list(nw_final=nw_final, baseline_final=baseline_final))
}






bestNW_dyn = function(n, obs, delta, nw_all, baselines_all, std, aucROC, aucROC_noNA, aucPR, t, k, graphSize_, outputDir, function_)
{
	# print best network for each value of std
	print("-------------------------------------------------------------")
	print(sprintf('t: %s, sd: %s outputDir: %s',t,std, outputDir))
	print(aucROC)
	
	print("obs")
	print(obs)
#	print(aucROC_noNA)
	
	# safeguard for when the program finds only empty nets
	if (length(aucROC_noNA) == 0){
		print("aucROC=noNA puff")
		nw_final = matrix(0, nrow=n, ncol=n)
		baseline_final = matrix(0, nrow=1, ncol=n)
		
		write.matrix(nw_final, file = sprintf("%s/adj_t=%s_k=%s_std=%s",outputDir,t,k,std))
		write.matrix(baseline_final, file = sprintf("%s/baseline_t=%s_k=%s_std=%s",outputDir,t,k,std))
		
		
		system(sprintf('python net5.py  printGraph_%s %s/adj_t=%s_k=%s_std=%s %s/k%s/std%s/graph%s_k%s_std%s',graphSize_, outputDir,t,k,std, outputDir,k,std,t,k,std), intern=TRUE)
		
		print(paste("max auc (zero): ", 0 , sep=""))
		print(nw_final)
		print(baseline_final)
	}
	else{
	
		ind_baseline = which(aucROC_noNA==max(aucROC_noNA,na.rm = TRUE))
		ind01 = n * (which(aucROC_noNA==max(aucROC_noNA,na.rm = TRUE))-1) + 1
		ind02 = ind01 + n -1
		nw_mean = baseline_mean= c()
		
		
		try(
			if (which(aucROC==max(aucROC,na.rm = TRUE)) != which(aucPR==max(aucPR,na.rm = TRUE))){
				print(" nw for best aucROC != nw for best aucPR")
			}, FALSE)
		
		
		for ( j in 1:length(ind01) ){
			nw_tmp = nw_all[ind01[j]:ind02[j],]
			nw_mean = rbind(nw_mean,c(t(nw_tmp)))
			baseline_tmp = baselines_all[ind_baseline[j],]
			baseline_mean = rbind(baseline_mean, baseline_tmp)
		}
		
	#	print(nw_all)
	#	print(nw_mean)
		nw_final = matrix(apply(nw_mean,2,median,na.rm=T),nrow=n,ncol=n,byrow=T)
		
		for (i in 1:n){
			for (j in 1:n){
				if (obs[j] >= delta[j]){
					nw_final[j,i] = obs[j] * nw_final[j,i]
				}
				else{
					nw_final[j,i] = 0
				}
			}
		}
		
		print("nw_final")
		print(nw_final)
		baseline_final = apply(baseline_mean, 2, median, na.rm=T)

		write.matrix(nw_final, file = sprintf("%s/adj_t=%s_k=%s_std=%s",outputDir,t,k,std))
		write.matrix(baseline_final, file = sprintf("%s/baseline_t=%s_k=%s_std=%s",outputDir,t,k,std))
		
		
		system(sprintf('python net5.py  printGraph_%s %s/adj_t=%s_k=%s_std=%s %s/k%s/std%s/graph%s_k%s_std%s',graphSize_, outputDir,t,k,std, outputDir,k,std,t,k,std), intern=TRUE)
		
		print(paste("max auc: ",max(aucROC_noNA,na.rm = TRUE), sep=""))
		print(nw_final)
		print(baseline_final)
	}
	return(list(nw_final=nw_final, baseline_final=baseline_final))
}





bestNWMSE = function(n, nw_all, T_nw, std, i, graphSize_, outputDir, function_)
{
	# print best network for each value of std
	print("-------------------------------------------------------------")
	print(sprintf('MSE t: %s, sd: %s outputDir: %s',i,std, outputDir))
#	print(aucROC)
#	print(aucROC_noNA)
	
#	print(nw_all)
	nw_final = getSampleMatrix(nw_all,n)
	
	colnames(nw_final) = NULL
#	print(nw_final)
	write.matrix(nw_final, file = sprintf("%s/MSEadj_t=%s_std=%s",outputDir,i, std))
	
	system(sprintf('python net%s.py  printGraph_%s %s/MSEadj_t=%s_std=%s %s/std%s/MSEgraph%s_%s', n,graphSize_, outputDir,i,std, outputDir,std,i,std), intern=TRUE)
	
	und_roc <- calc_ROC_underlying_nonIterative(c(t(nw_final)),T_nw,path="loocv_roc",triple=TRUE,plot=FALSE,t=t,std=std,outputDir=outputDir)

	print(paste("max auc roc ", und_roc[[2]], sep=""))
	print(paste("max auc pr ", und_roc[[3]], sep=""))
#	print(paste("max auc: ",max(aucROC_noNA,na.rm = TRUE), sep=""))
	print(nw_final)
#	print(baseline_final)
	
	return(list(nw_final=nw_final))
}






bestNWMSE_it = function(n, nw_all, T_nw, std, i, graphSize_, outputDir, function_)
{
	# print best network for each value of std
	print("-------------------------------------------------------------")
	print(sprintf('MSE t: %s, sd: %s outputDir: %s',i,std, outputDir))
#	print(aucROC)
#	print(aucROC_noNA)
	
#	print(nw_all)
	nw_final = getSampleMatrix(nw_all,n)
	
	colnames(nw_final) = NULL
#	print(nw_final)
	write.matrix(nw_final, file = sprintf("%s/MSEadj_t=%s_std=%s",outputDir,i, std))
	
	system(sprintf('python net%s.py  printGraph_%s %s/MSEadj_t=%s_std=%s %s/std%s/MSEgraph%s_%s', n,graphSize_, outputDir,i,std, outputDir,std,i,std), intern=TRUE)
	
	und_roc <- calc_ROC(nw_all,T_nw,path="loocv_roc",triple=TRUE,plot=FALSE,t=t,std=std,outputDir=outputDir)

	print(paste("max auc roc ", und_roc[[2]], sep=""))
	print(paste("max auc pr ", und_roc[[3]], sep=""))
#	print(paste("max auc: ",max(aucROC_noNA,na.rm = TRUE), sep=""))
	print(nw_final)
#	print(baseline_final)
	
	return(list(nw_final=nw_final))
}

