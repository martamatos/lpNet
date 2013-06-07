getBaseline = function(res,n,allpos=FALSE)
{

	if (allpos == FALSE) base=res$solution[(2*n*n+1):(2*n*n+n)]
	else base=res$solution[(n*n+1):(n*n+n)]
	
	return(base)
}
