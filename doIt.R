#
# returns the edges calculated using the best lambda
#
doIt <- function(function_,loocv_times,obs,n,b,K,delta,lambda=lamd,lambda2,annot,annot_node,T_=NULL,previousNet=NULL,baseline=NULL,previousBaseline=NULL,prior=NULL,startNode=NULL,endNode=NULL,allint=FALSE,active_mu,active_sd,inactive_mu,inactive_sd)
{

	lambda <- calcRangeLambda(delta,obs,stepsize=0.1)
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

	lambda <- calcRangeLambda(delta,obs,stepsize=0.1)
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



#
# returns the edges calculated using the best lambda
#
doIt_dyn <- function(function_,lambda2_,loocv_times,obs,n,b,K,delta,lambda=lamd,annot,annot_node,T_=NULL,previousNet=NULL,baseline=NULL,previousBaseline=NULL,prior=NULL,startNode=NULL,endNode=NULL,allint=FALSE,active_mu,active_sd,inactive_mu,inactive_sd)
{

	lambda <- calcRangeLambda_nonIterative(delta,obs,stepsize=0.05)
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
#
# returns the edges calculated using the best lambda
#
doIt_dynExp <- function(function_,lambda2_,loocv_times,obs,n,b,K,delta,lambda=lamd,annot,annot_node,T_=NULL,previousNet=NULL,baseline=NULL,previousBaseline=NULL,prior=NULL,startNode=NULL,endNode=NULL,allint=FALSE,active_mu,active_sd,inactive_mu,inactive_sd)
{

	lambda <- calcRangeLambda_nonIterative(delta,obs,stepsize=0.05)
	MSE <- Inf
  print("lambda")
  print(lambda)
	for(lamd in lambda)
	{
		loocv_res <- loocv_dyn(function_=function_,lambda2_=NULL,times=loocv_times,obs=obs,n=n,b=b,K=K,delta=delta,lambda=lamd,annot=annot,annot_node=annot_node,T_=T_,previousNet=NULL,baseline_=NULL,previousBaseline=NULL,active_mu=active_mu,active_sd=active_sd,inactive_mu=inactive_mu,inactive_sd=inactive_sd,prior=NULL,sourceNode=NULL,sinkNode=NULL,allint=FALSE,allpos=FALSE)
		
		if(loocv_res$MSE<MSE){
			MSE <- loocv_res$MSE
			edges_all <- loocv_res$edges_all
			baseline_all <- loocv_res$baseline_all
			bestLambda <- lamd
		}
	}
	print(paste("bestLambda ", bestLambda, sep=""))
	print(paste("bestMSE ", MSE, sep=""))
	print(edges_all)
	return(list(MSE=MSE,edges_all=edges_all,baseline_all=baseline_all,bestLambda=bestLambda))
}


doIt_exp <- function(function_,lambda2_,loocv_times,obs,n,b,K,delta,lambda=lamd,annot,annot_node,T_=NULL,previousNet=NULL,baseline=NULL,previousBaseline=NULL,prior=NULL,startNode=NULL,endNode=NULL,allint=FALSE,active_mu,active_sd,inactive_mu,inactive_sd)
{

	lambda <- calcRangeLambda(delta,obs,stepsize=0.05)
	MSE <- Inf
  print("lambda")
  print(lambda)
	for(lamd in lambda)
	{
		loocv_res <-loocv(function_=function_,times=loocv_times,obs=obs,n=n,b=b,K=K,delta=delta,lambda=lamd,lambda2=NULL,annot=getEdgeAnnot(n),annot_node=annot_node,T_=1,previousNet=NULL,baseline_=NULL,previousBaseline=NULL,active_mu=active_mu,active_sd=active_sd,inactive_mu=inactive_mu,inactive_sd=inactive_sd,prior=NULL,sourceNode=NULL,sinkNode=NULL,allint=FALSE,allpos=FALSE)
		
		if(loocv_res$MSE<MSE){
			MSE <- loocv_res$MSE
			edges_all <- loocv_res$edges_all
			baseline_all <- loocv_res$baseline_all
			bestLambda <- lamd
		}
	}
	print(paste("bestLambda ", bestLambda, sep=""))
	print(paste("bestMSE ", MSE, sep=""))
	print(edges_all)
	return(list(MSE=MSE,edges_all=edges_all,baseline_all=baseline_all,bestLambda=bestLambda))
}

doIt_dynExp_sahin03 <- function(function_,lambda2_,loocv_times,obs,n,b,K,delta,lambda=lamd,annot,annot_node,T_=NULL,previousNet=NULL,baseline=NULL,previousBaseline=NULL,prior=NULL,startNode=NULL,endNode=NULL,allint=FALSE,active_mu,active_sd,inactive_mu,inactive_sd)
{

	lambda <- calcRangeLambda_nonIterative_sahin03(delta,obs,stepsize=0.05)
	MSE <- Inf
  print("lambda")
  print(lambda)
	for(lamd in lambda)
	{
		loocv_res <- loocv_dyn_sahin03(function_=function_,lambda2_=NULL,times=loocv_times,obs=obs,n=n,b=b,K=K,delta=delta,lambda=lamd,annot=annot,annot_node=annot_node,T_=T_,previousNet=NULL,baseline_=NULL,previousBaseline=NULL,active_mu=active_mu,active_sd=active_sd,inactive_mu=inactive_mu,inactive_sd=inactive_sd,prior=NULL,sourceNode=NULL,sinkNode=NULL,allint=FALSE,allpos=FALSE)
		
		if(loocv_res$MSE<MSE){
			MSE <- loocv_res$MSE
			edges_all <- loocv_res$edges_all
			baseline_all <- loocv_res$baseline_all
			bestLambda <- lamd
		}
	}
	print(paste("bestLambda ", bestLambda, sep=""))
	print(paste("bestMSE ", MSE, sep=""))
	print(edges_all)
	return(list(MSE=MSE,edges_all=edges_all,baseline_all=baseline_all,bestLambda=bestLambda))
}

doIt_dynExp_ddepn <- function(function_,lambda2_,loocv_times,obs,n,b,K,delta,lambda=lamd,annot,annot_node,T_=NULL,previousNet=NULL,baseline=NULL,previousBaseline=NULL,prior=NULL,startNode=NULL,endNode=NULL,allint=FALSE,active_mu,active_sd,inactive_mu,inactive_sd)
{

	lambda <- calcRangeLambda_nonIterative_ddepn(delta,obs,stepsize=0.05)
	MSE <- Inf
  print("lambda")
  print(lambda)
	for(lamd in lambda)
	{
		loocv_res <- loocv_dyn_ddepn(function_=function_,lambda2_=NULL,times=loocv_times,obs=obs,n=n,b=b,K=K,delta=delta,lambda=lamd,annot=annot,annot_node=annot_node,T_=T_,previousNet=NULL,baseline_=NULL,previousBaseline=NULL,active_mu=active_mu,active_sd=active_sd,inactive_mu=inactive_mu,inactive_sd=inactive_sd,prior=NULL,sourceNode=NULL,sinkNode=NULL,allint=FALSE,allpos=FALSE)
		
		if(loocv_res$MSE<MSE){
			MSE <- loocv_res$MSE
			edges_all <- loocv_res$edges_all
			baseline_all <- loocv_res$baseline_all
			bestLambda <- lamd
		}
	}
	print(paste("bestLambda ", bestLambda, sep=""))
	print(paste("bestMSE ", MSE, sep=""))
	print(edges_all)
	return(list(MSE=MSE,edges_all=edges_all,baseline_all=baseline_all,bestLambda=bestLambda))
}


doIt_dynExp_ddepn_kfold <- function(function_,kfold,loocv_times,obs,n,b,K,delta,lambda=lamd,lambda2,annot,annot_node,T_=NULL,previousNet=NULL,baseline=NULL,previousBaseline=NULL,prior=NULL,startNode=NULL,endNode=NULL,allint=FALSE,active_mu,active_sd,inactive_mu,inactive_sd)
{

	lambda <- calcRangeLambda_nonIterative_ddepn(delta,obs,stepsize=0.1)
	MSE <- Inf
  print("lambda")
  print(lambda)
	for(lamd in lambda)
	{
		kfoldCV_res <- kfoldCV_dyn_ddepn(function_=function_,lambda2_=lambda2,kfold=kfold,times=loocv_times,obs=obs,n=n,b=b,K=K,delta=delta,lambda=lamd,annot=annot,annot_node=annot_node,T_=T_,previousNet=previousNet,baseline_=baseline_,previousBaseline=previousBaseline,active_mu=active_mu,active_sd=active_sd,inactive_mu=inactive_mu,inactive_sd=inactive_sd,prior=prior,sourceNode=startNode,sinkNode=endNode,allint=allint)
		
		if(kfoldCV_res$MSE<MSE){
			MSE <- kfoldCV_res$MSE
			edges_all <- kfoldCV_res$edges_all
			baseline_all <- kfoldCV_res$baseline_all
			bestLambda <- lamd
		}
	}
	print(paste("bestLambda ", bestLambda, sep=""))
	print(paste("bestMSE ", MSE, sep=""))
	print(edges_all)
	return(list(MSE=MSE,edges_all=edges_all,baseline_all=baseline_all,bestLambda=bestLambda))
}



doIt_dyn_kfold <- function(function_,kfold,loocv_times,obs,n,b,K,delta,lambda=lamd,lambda2,annot,annot_node,T_=NULL,previousNet=NULL,baseline=NULL,previousBaseline=NULL,prior=NULL,startNode=NULL,endNode=NULL,allint=FALSE,active_mu,active_sd,inactive_mu,inactive_sd)
{

	lambda <- calcRangeLambda_nonIterative(delta,obs,stepsize=0.1)
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

doIt_dyn_kfold005 <- function(function_,kfold,loocv_times,obs,n,b,K,delta,lambda=lamd,lambda2,annot,annot_node,T_=NULL,previousNet=NULL,baseline=NULL,previousBaseline=NULL,prior=NULL,startNode=NULL,endNode=NULL,allint=FALSE,active_mu,active_sd,inactive_mu,inactive_sd)
{

	lambda <- calcRangeLambda_nonIterative(delta,obs,stepsize=0.05)
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


doIt_dyn_kfold_2lev <- function(function_,kfold,loocv_times,obs,n,b,K,delta,lambda=lamd,lambda2,annot,annot_node,T_=NULL,previousNet=NULL,baseline=NULL,previousBaseline=NULL,prior=NULL,startNode=NULL,endNode=NULL,allint=FALSE,active_mu,active_sd,inactive_mu,inactive_sd,inactivePknock_mu,inactivePknock_sd)
{

	lambda <- calcRangeLambda_nonIterative(delta,obs,stepsize=0.05)
	MSE <- Inf
#  print("lambda")
#  print(lambda)
	for(lamd in lambda)
	{
		kfoldCV_res <- kfoldCV_dyn_2lev(function_=function_,lambda2_=lambda2,kfold=kfold,times=loocv_times,obs=obs,n=n,b=b,K=K,delta=delta,lambda=lamd,annot=annot,annot_node=annot_node,T_=T_,previousNet=previousNet,baseline_=baseline_,previousBaseline=previousBaseline,active_mu=active_mu,active_sd=active_sd,inactive_mu=inactive_mu,inactive_sd=inactive_sd,inactivePknock_mu=inactivePknock_mu, inactivePknock_sd=inactivePknock_sd, prior=prior,sourceNode=startNode,sinkNode=endNode,allint=allint)
		
		if(kfoldCV_res$MSE<MSE){
			MSE <- kfoldCV_res$MSE
			edges_all <- kfoldCV_res$edges_all
			baseline_all <- kfoldCV_res$baseline_all
			bestLambda <- lamd
		}
	}
	return(list(MSE=MSE,edges_all=edges_all,baseline_all=baseline_all,bestLambda=bestLambda))
}


doIt_dyn_kfoldExp <- function(function_,kfold,loocv_times,obs,n,b,K,delta,lambda=lamd,lambda2,annot,annot_node,T_=NULL,previousNet=NULL,baseline=NULL,previousBaseline=NULL,prior=NULL,startNode=NULL,endNode=NULL,allint=FALSE,active_mu,active_sd,inactive_mu,inactive_sd)
{

	lambda <- calcRangeLambda_nonIterative(delta,obs,stepsize=0.05)
	MSE <- Inf
  print("lambda")
  print(lambda)
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
	print(paste("bestLambda ", bestLambda, sep=""))
	print(paste("bestMSE ", MSE, sep=""))
	print(edges_all)
	return(list(MSE=MSE,edges_all=edges_all,baseline_all=baseline_all,bestLambda=bestLambda))
}

