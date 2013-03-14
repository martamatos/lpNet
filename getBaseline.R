getBaseline = function(res,n,T_=NULL,allpos=FALSE)
{
	if (length(T_)>0){
		if (allpos == FALSE) base=res$solution[(2*T_*n*n+1):(2*T_*n*n+n)]
		else base=res$solution[(T_*n*n+1):(T_*n*n+n)]
	}
	else{
		if (allpos == FALSE) base=res$solution[(2*n*n+1):(2*n*n+n)]
		else base=res$solution[(n*n+1):(n*n+n)]
	}
	
	return(base)
}
