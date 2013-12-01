options(digits=3)

source("../scripts/build_timeSeriesData.r")
source("../src/sourceDir.R")
sourceDir("../src")

test_bootstrap = function(){

	# setup problem

	#library(gurobi)
	library(lpSolve)
	
	LPfunction_lpSolve = get("doILP_dyn_discretized_dream8_new")
	LPfunction_gurobi = get("doILP_dyn_discretized_dream8_new_gurobi")
	getAdja_function= get("getAdja")
	getBaseline_function= get("getBaseline")
#	getAdja_function= get("getAdja_gurobi")
#	getBaseline_function= get("getBaseline_gurobi")
	
	n = 3 # number of genes
	K = 2 # number of knockdowns

	annot_node <- seq(1,n)
	annot <- getEdgeAnnot(n)

	## calculate Observation matrix with different Gaussian for activation and deactivation
	active_mu <- 0.95
	inactive_mu <- 0.56
	active_sd <- 0.01
	inactive_sd <- 0.01
	# use three replicates
	replnum <- 3
	
	bootstrap_times = 10
	
	sd_all <- c(0.05) # use different standard deviations for different noise levels

	kfold=  NULL

	muPgene = FALSE
	muPgk = FALSE
	muPgt = FALSE
	muPgkt = FALSE

	deltaPk = FALSE
	deltaPt = FALSE
	deltaPkt = FALSE

	#-------------------------------------
	# define network data
	#-------------------------------------

	T_undNet = matrix(c(0,1,0,
											0,0,1,
											0,0,0), nrow=n,ncol=n,byrow=T)
										
	nodeNames = seq(1,n)
	
	bvec = c(1,1,1,
					 1,0,1)
	
	timeSeriesData = build_timeSeries(T_undNet,nodeNames,bvec,n,K)
	print(timeSeriesData)
	T_ = timeSeriesData$T_
	geneState =timeSeriesData$geneStateVec
	
	act_mat = calcActivation(T_undNet,bvec,n,K)
	print(act_mat)
	
	print(geneState)
	sd_i = 1

	# generate observation matrix and delta vector
	obs = array(NA, c(n,K,T_, replnum))
	tmp = c()
	
	for (t in 1:T_){
		geneState_ = geneState[,,t]
		obs_all <-  vector()

		for(repl in 1:replnum)
		{
			obs_tmp <- getObsMat(act_mat,active_mu,sd_all[1],inactive_mu,sd_all[1],geneState_)
			obs[,,t,repl] = obs_tmp
			obs_all <- rbind(obs_all,c(obs_tmp))
			tmp <- cbind(tmp,rnorm(n,mean(c(active_mu,inactive_mu)),sd_all[1]))
		}
	}
	delta <- apply(tmp,1,mean,na.rm=T)

#	obs = array(NA, c(n,K,T_))
#	for (t in 1:T_){
#		geneState_ = geneState[,,t]
#		obs_all <-  vector()

#		for(repl in 1:replnum)
#		{
#			obs_tmp <- getObsMat(act_mat,active_mu,sd_all[1],inactive_mu,sd_all[1],geneState_)
#			obs_all <- rbind(obs_all,c(obs_tmp))
#			tmp <- cbind(tmp,rnorm(n,mean(c(active_mu,inactive_mu)),sd_all[1]))
#		}
#		obs[,,t] <- matrix(apply(obs_all,2,mean,na.rm=T),nrow=n,ncol=K)
#	}
#	delta <- apply(tmp,1,mean,na.rm=T)

	print(obs)


	res = bootstrap_dyn(function_=LPfunction_lpSolve,getAdja_function, getBaseline_function,predFunction=NULL,kfold=NULL,times=bootstrap_times,obs=obs,n=n,b=bvec,K=K,delta=delta,lambda=0.1,annot=annot,annot_node=annot_node,T_=T_,active_mu=active_mu,active_sd=active_sd,inactive_mu=inactive_mu,inactive_sd=inactive_sd,prior=NULL,sourceNode=NULL,sinkNode=NULL,allint=FALSE,allpos=FALSE,muPgene=FALSE,muPgk=FALSE,muPgt=FALSE,muPgkt=FALSE,deltaPk=FALSE,deltaPt=FALSE,deltaPkt=FALSE)
	
#	print(res)
#	loocv_dyn(function_=LPfunction_lpSolve,predFunction=get("calcPredictionLOOCV_dyn_disc"),getAdja_function=getAdja_function, getBaseline_function=getBaseline_function,kfold=NULL,times=bootstrap_times,obs=obs,n=n,b=bvec,K=K,delta=delta,lambda=0.1,annot=annot,annot_node=annot_node,T_=T_,active_mu=active_mu,active_sd=sd_all[1],inactive_mu=inactive_mu,inactive_sd=sd_all[1],prior=NULL,sourceNode=NULL,sinkNode=NULL,allint=FALSE,allpos=FALSE,muPgene=FALSE,muPgk=FALSE,muPgt=FALSE,muPgkt=FALSE,deltaPk=FALSE,deltaPt=FALSE,deltaPkt=FALSE)
	
#	(function_,predFunction=NULL,getAdja_function, getBaseline_function,kfold=NULL,times,obs,n,b,K,delta,lambda,annot,annot_node,T_,active_mu=NULL,active_sd=NULL,inactive_mu=NULL,inactive_sd=NULL,prior=NULL,sourceNode=NULL,sinkNode=NULL,allint=FALSE,allpos=FALSE,muPgene=FALSE,muPgk=FALSE,muPgt=FALSE,muPgkt=FALSE,deltaPk=FALSE,deltaPt=FALSE,deltaPkt=FALSE)

}
test_bootstrap()
