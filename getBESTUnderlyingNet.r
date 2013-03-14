getBESTUnderlyingNet = function(n, nw_all, std, aucROC_und, aucPR_und, graphSize_, outputDir, function_)
{

	# print best underlying network for each value of std, calculaed from all the underlying networks for each step
	print("-------------------------------------------------------------")
	print(sprintf('Best underlying network each step, sd: %s outputDir: %s',std, outputDir))
#	print(nw_all)
#	print(aucROC_und)
#	print(aucPR_und)
	ind01 = which(aucROC_und==max(aucROC_und,na.rm = TRUE))
	nw_mean = c()
	
	try(
		if (which(aucROC_und==max(aucROC_und,na.rm = TRUE)) != which(aucPR_und==max(aucPR_und,na.rm = TRUE))){
			print(" nw for best aucROC != nw for best aucPR")
		}, FALSE)
	

	for ( j in 1:length(ind01) ){
		nw_tmp = nw_all[ind01[j],]
		nw_mean = rbind(nw_mean,c(t(nw_tmp)))
	}
	
#	print(nw_mean)
	nw_final = matrix(apply(nw_mean,2,median,na.rm=T),nrow=n,ncol=n,byrow=T)
	
	write.matrix(nw_final, file = sprintf("%s/underStep_std=%s",outputDir,std))
	
	system(sprintf('python net%s.py  printGraph_%s %s/underStep_std=%s %s/graph_underStep_%s', n,graphSize_, outputDir, std, outputDir, std), intern=TRUE)
	
	print(paste("max auc: ",max(aucROC_und,na.rm = TRUE), sep=""))
	print(nw_final)
	
	return(nw_final)
}


getBESTUnderlyingNetMSE = function(n, nw_all, std, aucROC_und, aucPR_und, graphSize_, outputDir, function_)
{

	# print best underlying network for each value of std, calculaed from all the underlying networks for each step
	print("-------------------------------------------------------------")
	print(sprintf('Best underlying network each step MSE, sd: %s outputDir: %s',std, outputDir))
#	print(nw_all)
#	print(aucROC_und)
#	print(aucPR_und)
	ind01 = which(aucROC_und==max(aucROC_und,na.rm = TRUE))
	nw_mean = c()
	
	
#	print("lenght")
#	print(length(aucROC_und))
#	print(length(nw_all))
#	print(dim(nw_all))
	
	try(
		if (which(aucROC_und==max(aucROC_und,na.rm = TRUE)) != which(aucPR_und==max(aucPR_und,na.rm = TRUE))){
			print(" nw for best aucROC != nw for best aucPR")
		}, FALSE)
	

	for ( j in 1:length(ind01) ){
		nw_tmp = nw_all[ind01[j],]
		nw_mean = rbind(nw_mean,c(t(nw_tmp)))
	}
#	print(nw_mean)
	nw_final = matrix(apply(nw_mean,2,median,na.rm=T),nrow=n,ncol=n, byrow=T)
#	print(nw_final)
	write.matrix(nw_final, file = sprintf("%s/underStepMSE_std=%s",outputDir,std))
	
	system(sprintf('python net%s.py  printGraph_%s %s/underStepMSE_std=%s %s/graph_underStepMSE_%s', n,graphSize_, outputDir, std, outputDir, std), intern=TRUE)
	
	print(paste("max auc: ",max(aucROC_und,na.rm = TRUE), sep=""))
	print(nw_final)
	
	return(nw_final)
}

