#
# returns the edges calculated using the best lambda
#
doIt_LP <- function(LPfunction,CVfunction, predFunction, loocv_times, kfold=NULL,stepsize,obs,n,b,K,delta,annot,annot_node,T_=NULL,prior=NULL,startNode=NULL,endNode=NULL,allint=FALSE,allpos=FALSE,active_mu,active_sd,inactive_mu,inactive_sd,muPgene=FALSE,muPgk=FALSE,muPgt=FALSE,muPgkt=FALSE,deltaPk=FALSE,deltaPt=FALSE,deltaPkt=FALSE)
{

	lambda <- calcRangeLambda_LP(delta,obs,stepsize=stepsize,deltaPk=deltaPk)
	MSE <- Inf
	
  print("all lambdas")
  print(lambda)
  
	for(lamd in lambda)
	{
		print(paste("current lambda: ", lamd, sep=""))
		
		res <- CVfunction(function_=LPfunction,predFunction=predFunction,kfold=kfold,times=loocv_times,obs=obs,n=n,b=b,K=K,delta=delta,lambda=lamd,annot=annot,annot_node=annot_node,T_=T_,active_mu=active_mu,active_sd=active_sd,inactive_mu=inactive_mu,inactive_sd=inactive_sd,prior=prior,sourceNode=startNode,sinkNode=endNode,allint=allint,allpos=allpos,muPgene=muPgene,muPgk=muPgk,muPgt=muPgt,muPgkt=muPgkt,deltaPk=deltaPk,deltaPt=deltaPt,deltaPkt=deltaPkt)
		
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
doIt_dyn<- function(LPfunction,CVfunction, predFunction, loocv_times, kfold=NULL,stepsize,obs,n,b,K,delta,annot,annot_node,T_=NULL,prior=NULL,startNode=NULL,endNode=NULL,allint=FALSE,allpos=FALSE,active_mu,active_sd,inactive_mu,inactive_sd,muPgene=FALSE,muPgk=FALSE,muPgt=FALSE,muPgkt=FALSE,deltaPk=FALSE,deltaPt=FALSE,deltaPkt=FALSE)
{

	lambda <- calcRangeLambda_dyn(delta,obs,stepsize=stepsize,deltaPk=deltaPk,deltaPt=deltaPt,deltaPkt=deltaPkt)
	MSE <- Inf
	
  print("all lambdas")
  print(lambda)
  
	for(lamd in lambda)
	{
		print(paste("current lambda: ", lamd, sep=""))
		
		res <- CVfunction(function_=LPfunction,predFunction=predFunction,kfold=kfold,times=loocv_times,obs=obs,n=n,b=b,K=K,delta=delta,lambda=lamd,annot=annot,annot_node=annot_node,T_=T_,active_mu=active_mu,active_sd=active_sd,inactive_mu=inactive_mu,inactive_sd=inactive_sd,prior=prior,sourceNode=startNode,sinkNode=endNode,allint=allint,allpos=allpos,muPgene=muPgene,muPgk=muPgk,muPgt=muPgt,muPgkt=muPgkt,deltaPk=deltaPk,deltaPt=deltaPt,deltaPkt=deltaPkt)

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












# not in use


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


