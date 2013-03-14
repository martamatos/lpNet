getUnderlyingNet_indiv = function(net, T_, n, std, outputDir, graphSize_)
{
#	print("getUnderlyingNet_indiv")
#	print(net)
	tempUndNet = net
	
	tempUndNet_median = apply(tempUndNet, 2, median)
	tempUndNet_mad = apply(tempUndNet, 2, mad)
	
	
	underlyingNet = c(rep(0,n*n))
	underlyingNet_ind = which(abs(tempUndNet_median)>abs(tempUndNet_mad), arr.ind=T)
	underlyingNet[underlyingNet_ind] = tempUndNet_median[underlyingNet_ind]
#	underlyingNet = matrix(c(underlyingNet), nrow = n, ncol=n, byrow=T)
	
#	print("und net indiv")
#	print(tempUndNet_median)
#	print(tempUndNet_mad)
#	print(underlyingNet)
	return(underlyingNet)
}


getUnderlyingNet_indivMSE = function(net, T_, n, std, outputDir, graphSize_)
{
#	print(" net MSE")
#	print(net)
	tempUndNet = c()
	for (t in 1:T_){
		tempUndNet = rbind(tempUndNet, c(t(getSampleMatrix(net[[t]],n))))
	}
#	print("temp und net MSE")
#	print(dim(tempUndNet))
#	print(tempUndNet)
	
	
	tempUndNet_median = apply(tempUndNet, 2, median)
	tempUndNet_mad = apply(tempUndNet, 2, mad)
	
	

	underlyingNet = c(rep(0,n*n))
	underlyingNet_ind = which(abs(tempUndNet_median)>abs(tempUndNet_mad), arr.ind=T)
	underlyingNet[underlyingNet_ind] = tempUndNet_median[underlyingNet_ind]
#	underlyingNet = matrix(c(underlyingNet), nrow = n, ncol=n, byrow=T)
	
#	print("und net indivMSE")
#	print(tempUndNet_median)
#	print(tempUndNet_mad)
#	print(underlyingNet)
	return(underlyingNet)
}
