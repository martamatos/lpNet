getUnderlyingNet = function(active_networks, T_undNet, T_, n, std, outputDir, graphSize_)
{

	tempUndNet = c()
#	print("getUnderlyingNet")
#	print(active_networks)
	for (t in 1:T_){
		tempUndNet = rbind(tempUndNet, c(t(active_networks[[t]]$nw_final)))
	}
#	print(active_networks)
#	print(tempUndNet)
	
	tempUndNet_median = apply(tempUndNet, 2, median)
	tempUndNet_mad = apply(tempUndNet, 2, mad)
	
#	for (i in 1:(n*n)){
#		print(tempUndNet_median[i])
#		print(tempUndNet_mad[i])
#		print("------")
#	}

	underlyingNet = c(rep(0,n*n))
	underlyingNet_ind = which(abs(tempUndNet_median)>abs(tempUndNet_mad), arr.ind=T)
	underlyingNet[underlyingNet_ind] = tempUndNet_median[underlyingNet_ind]
	
	underlyingNet2 = matrix(c(underlyingNet), nrow = n, ncol=n, byrow=T)
#	print(underlyingNet2)
	write.matrix(underlyingNet2, file = sprintf("%s/underNet_std=%s",outputDir, std))
	system(sprintf('python net%s.py  printGraph_%s %s/underNet_std=%s %s/graphUnderNet_%s', n,graphSize_, outputDir,std, outputDir,std), intern=TRUE)
	
	und_roc <- calc_ROC_underlying_nonIterative(underlyingNet,T_undNet,path="loocv_roc",triple=TRUE,plot=FALSE,t=t,std=std,outputDir=outputDir)


	print("-------------------------------------------------------------")
	print(sprintf('Best underlying network efrom best active networks, sd: %s outputDir: %s',std, outputDir))
	print(underlyingNet2)
	print(paste("max auc roc ", und_roc[[2]], sep=""))
	print(paste("max auc pr ", und_roc[[3]], sep=""))
	
#	return(underlyingNet)
}



getUnderlyingNetMSE = function(active_networks, T_undNet, T_, n, std, outputDir, graphSize_)
{
	tempUndNet = c()
#	print("getUnderlyingNetMSE")
#	print(dim(active_networks))
#	print(active_networks)
	
	for (t in 1:T_){
		tempUndNet = rbind(tempUndNet, c(t(getSampleMatrix(active_networks[[t]],n))))
	}
#	print(active_networks)
#	print(tempUndNet)
	
	tempUndNet_median = apply(tempUndNet, 2, median)
	tempUndNet_mad = apply(tempUndNet, 2, mad)
	
#	for (i in 1:(n*n)){
#		print(tempUndNet_median[i])
#		print(tempUndNet_mad[i])
#		print("------")
#	}

	underlyingNet = c(rep(0,n*n))
	underlyingNet_ind = which(abs(tempUndNet_median)>abs(tempUndNet_mad), arr.ind=T)
	underlyingNet[underlyingNet_ind] = tempUndNet_median[underlyingNet_ind]

	underlyingNet2 = matrix(c(underlyingNet), nrow=n, ncol=n, byrow=T)
#	print("length(underlyingNet)")
#	print(length(underlyingNet2))
#	print(underlyingNet2)
	
	write.matrix(underlyingNet2, file = sprintf("%s/underNetMSE_std=%s",outputDir, std))
	system(sprintf('python net%s.py  printGraph_%s %s/underNetMSE_std=%s %s/graphUnderNetMSE_%s', n,graphSize_, outputDir,std, outputDir,std), intern=TRUE)
	
	und_roc <- calc_ROC_underlying_nonIterative(underlyingNet,T_undNet,path="loocv_roc",triple=TRUE,plot=FALSE,t=t,std=std,outputDir=outputDir)

	print("-------------------------------------------------------------")
	print(sprintf('Best underlying network efrom best active networks MSE, sd: %s outputDir: %s',std, outputDir))
	print(underlyingNet2)
	print(paste("max auc roc ", und_roc[[2]], sep=""))
	print(paste("max auc pr ", und_roc[[3]], sep=""))
	
#	return(underlyingNet)
}
