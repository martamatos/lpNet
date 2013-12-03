#
# returns the edges calculated using the best lambda
#
doIt_LP <- function(LPfunction,CVfunction, predFunction,getAdja_function, getBaseline_function, loocv_times, kfold=NULL,stepsize,obs,n,b,K,delta,annot,annot_node,T_=NULL,prior=NULL,startNode=NULL,endNode=NULL,allint=FALSE,allpos=FALSE,active_mu,active_sd,inactive_mu,inactive_sd,mu_type=mu_type,delta_type=delta_type)
{

	lambda <- calcRangeLambda_LP(delta,obs,stepsize=stepsize,delta_type=delta_type)
	MSE <- Inf
	
  print("all lambdas")
  print(lambda)
  
	for(lamd in lambda)
	{
		print(paste("current lambda: ", lamd, sep=""))
		
		res <- CVfunction(function_=LPfunction,predFunction=predFunction,getAdja_function=getAdja_function, getBaseline_function=getBaseline_function,kfold=kfold,times=loocv_times,obs=obs,n=n,b=b,K=K,delta=delta,lambda=lamd,annot=annot,annot_node=annot_node,T_=T_,active_mu=active_mu,active_sd=active_sd,inactive_mu=inactive_mu,inactive_sd=inactive_sd,prior=prior,sourceNode=startNode,sinkNode=endNode,allint=allint,allpos=allpos,mu_type=mu_type,delta_type=delta_type)
		
		if(res$MSE<MSE){
			MSE <- res$MSE
			edges_all <- res$edges_all
			baseline_all <- res$baseline_all
			bestLambda <- lamd
		}
	}
	
	print(paste("bestLambda ", bestLambda, sep=""))
	print(paste("bestMSE ", MSE, sep=""))
	
	return(list(MSE=MSE,edges_all=edges_all,baseline_all=baseline_all,bestLambda=bestLambda))
}

#
# returns the edges calculated using the best lambda
#
doIt_dyn<- function(LPfunction,CVfunction, predFunction,getAdja_function, getBaseline_function, realTime_function=NULL, loocv_times, kfold=NULL,stepsize,obs,n,b,K,delta,annot,annot_node,T_=NULL,prior=NULL,startNode=NULL,endNode=NULL,allint=FALSE,allpos=FALSE,active_mu,active_sd,inactive_mu,inactive_sd,mu_type=mu_type,delta_type=delta_type)
{

	lambda <- calcRangeLambda_dyn(delta,obs,stepsize=stepsize,delta_type=delta_type)
	MSE <- Inf
	
  print("all lambdas")
  print(lambda)
  
	for(lamd in lambda)
	{
		print(paste("current lambda: ", lamd, sep=""))
		
		res <- CVfunction(function_=LPfunction,predFunction=predFunction,getAdja_function=getAdja_function, getBaseline_function=getBaseline_function,realTime_function=realTime_function,kfold=kfold,times=loocv_times,obs=obs,n=n,b=b,K=K,delta=delta,lambda=lamd,annot=annot,annot_node=annot_node,T_=T_,active_mu=active_mu,active_sd=active_sd,inactive_mu=inactive_mu,inactive_sd=inactive_sd,prior=prior,sourceNode=startNode,sinkNode=endNode,allint=allint,allpos=allpos,mu_type=mu_type,delta_type=delta_type)

		if(res$MSE<MSE){
			MSE <- res$MSE
			edges_all <- res$edges_all
			baseline_all <- res$baseline_all
			bestLambda <- lamd
		}
	}
	
	print(paste("bestLambda ", bestLambda, sep=""))
	print(paste("bestMSE ", MSE, sep=""))
	
	return(list(MSE=MSE,edges_all=edges_all,baseline_all=baseline_all,bestLambda=bestLambda))
}



doIt_dyn_lambdaRange<- function(LPfunction,CVfunction, predFunction,getAdja_function, getBaseline_function,realTime_function=NULL, loocv_times, kfold=NULL,stepsize,obs,n,b,K,delta,annot,annot_node,T_=NULL,prior=NULL,startNode=NULL,endNode=NULL,allint=FALSE,allpos=FALSE,active_mu,active_sd,inactive_mu,inactive_sd,mu_type=mu_type,delta_type=delta_type)
{

	lambda <- calcRangeLambda_dyn_dream8(delta,obs,stepsize=stepsize,delta_type=delta_type)
	MSE <- Inf
	
  print("all lambdas")
  print(lambda)
  
	for(lamd in lambda)
	{
		print(paste("current lambda: ", lamd, sep=""))
		
		res <- CVfunction(function_=LPfunction,predFunction=predFunction,getAdja_function=getAdja_function, getBaseline_function=getBaseline_function,realTime_function=realTime_function,kfold=kfold,times=loocv_times,obs=obs,n=n,b=b,K=K,delta=delta,lambda=lamd,annot=annot,annot_node=annot_node,T_=T_,active_mu=active_mu,active_sd=active_sd,inactive_mu=inactive_mu,inactive_sd=inactive_sd,prior=prior,sourceNode=startNode,sinkNode=endNode,allint=allint,allpos=allpos,mu_type=mu_type,delta_type=delta_type)

		if(res$MSE<MSE){
			MSE <- res$MSE
			edges_all <- res$edges_all
			baseline_all <- res$baseline_all
			bestLambda <- lamd
		}
	}
	
	print(paste("bestLambda ", bestLambda, sep=""))
	print(paste("bestMSE ", MSE, sep=""))
	
	return(list(MSE=MSE,edges_all=edges_all,baseline_all=baseline_all,bestLambda=bestLambda))
}





