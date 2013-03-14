#
# returns the edges calculated using the best lambda
#
doIt <- function(function_,loocv_times,obs,n,b,K,delta,lambda=lamd,lambda2,annot,annot_node,T_=NULL,previousNet=NULL,baseline=NULL,previousBaseline=NULL,prior=NULL,startNode=NULL,endNode=NULL,allint=FALSE,active_mu,active_sd,inactive_mu,inactive_sd)
{

	lambda <- calcRangeLambda(delta,obs,stepsize=0.01)
	MSE <- Inf
#  print("lambda")
#  print(lambda)
	for(lamd in lambda)
	{
		loocv_res <- loocv(function_=function_,times=loocv_times,obs=obs,n=n,b=b,K=K,delta=delta,lambda=lamd,lambda2=lambda2,annot=getEdgeAnnot(n),annot_node=annot_node,T_=T_,previousNet=previousNet,baseline_=baseline,previousBaseline=previousBaseline,active_mu=active_mu,active_sd=active_sd,inactive_mu=inactive_mu,inactive_sd=inactive_sd,prior=NULL,sourceNode=NULL,sinkNode=NULL,allint=FALSE,allpos=FALSE)
			
		if(loocv_res$MSE<MSE){
			MSE <- loocv_res$MSE
			edges_all <- loocv_res$edges_all
			baseline_all <- loocv_res$baseline_all
			bestLambda <- lamd
		}
	}
	
	print(paste("bestLambda ", bestLambda, sep=""))
	print(paste("bestMSE ", MSE, sep=""))
	
	return(list(MSE=MSE,edges_all=edges_all,baseline_all=baseline_all,bestLambda=bestLambda))
}

doIt_kfold <- function(function_,kfold,loocv_times,obs,n,b,K,delta,lambda=lamd,lambda2,annot,annot_node,T_=NULL,previousNet=NULL,baseline=NULL,previousBaseline=NULL,prior=NULL,startNode=NULL,endNode=NULL,allint=FALSE,active_mu,active_sd,inactive_mu,inactive_sd)
{

	lambda <- calcRangeLambda(delta,obs,stepsize=0.01)
	MSE <- Inf
#  print("lambda")
#  print(lambda)
	for(lamd in lambda)
	{
		kfoldCV_res <- kfoldCV(function_=function_,lambda2_=lambda2,kfold=kfold,times=loocv_times,obs=obs,n=n,b=b,K=K,delta=delta,lambda=lamd,annot=annot,annot_node=annot_node,T_=T_,previousNet=previousNet,baseline_=baseline_,previousBaseline=previousBaseline,active_mu=active_mu,active_sd=active_sd,inactive_mu=inactive_mu,inactive_sd=inactive_sd,prior=prior,sourceNode=startNode,sinkNode=endNode,allint=allint)
		
		if(kfoldCV_res$MSE<MSE){
			MSE <- kfoldCV_res$MSE
			edges_all <- kfoldCV_res$edges_all
			baseline_all <- kfoldCV_res$baseline_all
			bestLambda <- lamd
		}
	}
	return(list(MSE=MSE,edges_all=edges_all,baseline_all=baseline_all,bestLambda=bestLambda))
}


doIt_lambda2 <- function(function_,lambda2_,loocv_times,obs,n,b,K,delta,lambda=lamd,annot,annot_node,T_=NULL,previousNet=NULL,baseline=NULL,previousBaseline=NULL,prior=NULL,startNode=NULL,endNode=NULL,allint=FALSE,active_mu,active_sd,inactive_mu,inactive_sd)
{

	lambda <- calcRangeLambda(delta,obs,stepsize=0.05)
	lambda2List = c(0.0, 0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2., 3., 5., 10., 100.)
	MSE <- Inf
#  print("lambda")
#  print(lambda)
	for(lamd in lambda)
	{
		for (lambda2_ in lambda2List)
		{
			loocv_res <- loocv(function_=function_,times=loocv_times,obs=obs,n=n,b=b,K=K,delta=delta,lambda=lamd,lambda2=lambda2_,annot=getEdgeAnnot(n),annot_node=annot_node,T_=T_,previousNet=previousNet,baseline_=baseline,previousBaseline=previousBaseline,active_mu=active_mu,active_sd=active_sd,inactive_mu=inactive_mu,inactive_sd=inactive_sd,prior=NULL,sourceNode=NULL,sinkNode=NULL,allint=FALSE,allpos=FALSE)
		
			if(loocv_res$MSE<MSE){
				MSE <- loocv_res$MSE
				edges_all <- loocv_res$edges_all
				baseline_all <- loocv_res$baseline_all
				bestLambda <- lamd
				bestLambda2 <- lambda2_
			}
		}
	}
	print(paste("bestLambda ", bestLambda, sep=""))
	print(paste("bestLambda2 ", bestLambda2, sep=""))
	print(paste("bestMSE ", MSE, sep=""))
	return(list(MSE=MSE,edges_all=edges_all,baseline_all=baseline_all,bestLambda=bestLambda))
}




doIt_lambda2_kfold <- function(function_,kfold,loocv_times,obs,n,b,K,delta,lambda=lamd,lambda2,annot,annot_node,T_=NULL,previousNet=NULL,baseline=NULL,previousBaseline=NULL,prior=NULL,startNode=NULL,endNode=NULL,allint=FALSE,active_mu,active_sd,inactive_mu,inactive_sd)
{

	lambda <- calcRangeLambda(delta,obs,stepsize=0.05)
	lambda2List = c(0.0, 0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2., 3., 5., 10., 100.)
	MSE <- Inf
#  print("lambda")
#  print(lambda)
	for(lamd in lambda)
	{
		for (lambda2_ in lambda2List)
		{
			kfoldCV_res <- kfoldCV(function_=function_,lambda2_=lambda2_,kfold=kfold,times=loocv_times,obs=obs,n=n,b=b,K=K,delta=delta,lambda=lamd,annot=annot,annot_node=annot_node,T_=T_,previousNet=previousNet,baseline_=baseline_,previousBaseline=previousBaseline,active_mu=active_mu,active_sd=active_sd,inactive_mu=inactive_mu,inactive_sd=inactive_sd,prior=prior,sourceNode=startNode,sinkNode=endNode,allint=allint)
			
			if(kfoldCV_res$MSE<MSE){
				MSE <- kfoldCV_res$MSE
				edges_all <- kfoldCV_res$edges_all
				baseline_all <- kfoldCV_res$baseline_all
				bestLambda <- lamd
				bestLambda2 <- lambda2_
			}
		}
	}
	print(paste("bestLambda ", bestLambda, sep=""))
	print(paste("bestLambda2 ", bestLambda2, sep=""))
	print(paste("bestMSE ", MSE, sep=""))
	return(list(MSE=MSE,edges_all=edges_all,baseline_all=baseline_all,bestLambda=bestLambda))
}


#
# returns the edges calculated using the best lambda
#
doIt_dyn <- function(function_,lambda2_,loocv_times,obs,n,b,K,delta,lambda=lamd,annot,annot_node,T_=NULL,previousNet=NULL,baseline=NULL,previousBaseline=NULL,prior=NULL,startNode=NULL,endNode=NULL,allint=FALSE,active_mu,active_sd,inactive_mu,inactive_sd)
{

	lambda <- calcRangeLambda_nonIterative(delta,obs,stepsize=0.01)
	MSE <- Inf
#  print("lambda")
#  print(lambda)
	for(lamd in lambda)
	{
		loocv_res <- loocv_dyn(function_,lambda2_=NULL,loocv_times,obs,n,b,K,delta,lambda=lamd,annot=getEdgeAnnot(n),annot_node,T_,previousNet,baseline,previousBaseline,active_mu,active_sd,inactive_mu,inactive_sd)
		
		if(loocv_res$MSE<MSE){
			MSE <- loocv_res$MSE
			edges_all <- loocv_res$edges_all
			baseline_all <- loocv_res$baseline_all
			bestLambda <- lamd
		}
	}
	print(paste("bestLambda ", bestLambda, sep=""))
	print(paste("bestMSE ", MSE, sep=""))
	return(list(MSE=MSE,edges_all=edges_all,baseline_all=baseline_all,bestLambda=bestLambda))
}



doIt_dyn_kfold <- function(function_,kfold,loocv_times,obs,n,b,K,delta,lambda=lamd,lambda2,annot,annot_node,T_=NULL,previousNet=NULL,baseline=NULL,previousBaseline=NULL,prior=NULL,startNode=NULL,endNode=NULL,allint=FALSE,active_mu,active_sd,inactive_mu,inactive_sd)
{

	lambda <- calcRangeLambda_nonIterative(delta,obs,stepsize=0.01)
	MSE <- Inf
#  print("lambda")
#  print(lambda)
	for(lamd in lambda)
	{
		kfoldCV_res <- kfoldCV_dyn(function_=function_,lambda2_=lambda2,kfold=kfold,times=loocv_times,obs=obs,n=n,b=b,K=K,delta=delta,lambda=lamd,annot=annot,annot_node=annot_node,T_=T_,previousNet=previousNet,baseline_=baseline_,previousBaseline=previousBaseline,active_mu=active_mu,active_sd=active_sd,inactive_mu=inactive_mu,inactive_sd=inactive_sd,prior=prior,sourceNode=startNode,sinkNode=endNode,allint=allint)
		
		if(kfoldCV_res$MSE<MSE){
			MSE <- kfoldCV_res$MSE
			edges_all <- kfoldCV_res$edges_all
			baseline_all <- kfoldCV_res$baseline_all
			bestLambda <- lamd
		}
	}
	return(list(MSE=MSE,edges_all=edges_all,baseline_all=baseline_all,bestLambda=bestLambda))
}

