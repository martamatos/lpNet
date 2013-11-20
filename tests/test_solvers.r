.setUp = function(){
	a = 2
}


test_gurobi_agains_lpSolve = function(){

	# setup problem
	
	library(gurobi)
	library(lpSolve)
	
	LPfunction_lpSolve = get("doILP_dyn_discretized_dream8_new")
	LPfunction_gurobi = get("doILP_dyn_discretized_dream8_new_gurobi")

	n = 16 # number of genes
	K = 17 # number of knockdowns

	annot_node <- seq(1,n)
	annot <- getEdgeAnnot(n)

	## calculate Observation matrix with different Gaussian for activation and deactivation
	active_mu <- 0.95
	inactive_mu <- 0.56
	active_sd <- 0.01
	inactive_sd <- 0.01
	# use three replicates
	replnum <- 3

	sd_all <- c(0.05) # use different standard deviations for different noise levels
	
	stepsize = 0.1
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

	T_undNet = matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
											0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
											0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,
											0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,1,
											0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
											1,1,0,0,0,0,0,0,0,0,1,0,0,1,0,0,
											0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,
											0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
											0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
											0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
											0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
											0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
											0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
											0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
											0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
											0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), nrow=n,ncol=n,byrow=T)
										
	nodeNames = seq(1,n)
	bvec = rep(1, n*K)
	for (i in 1:n){
		bvec[i*n+i] = 0
	}

	timeSeriesData = build_timeSeriesData(T_undNet,nodeNames,bvec,n,K)
	T_ = timeSeriesData$T_
	geneState =timeSeriesData$geneStateVec
	act_mat = calcActivation_dyn(T_undNet,bvec,n,K)
	
	sd_i = 1


	# generate observation matrix and delta vector
	obs = array(NA, c(n,K,T_))
	tmp = c()
	
	for (t in 1:T_){
		geneState_ = geneState[,,t]
		obs_all <-  vector()

		for(repl in 1:replnum)
		{
			obs_tmp <- getObsMat(act_mat,active_mu,sd_all[1],inactive_mu,sd_all[1],geneState_)
			obs_all <- rbind(obs_all,c(obs_tmp))
			tmp <- cbind(tmp,rnorm(n,mean(c(active_mu,inactive_mu)),sd_all[1]))
		}
		obs[,,t] <- matrix(apply(obs_all,2,mean,na.rm=T),nrow=n,ncol=K)
	}
	delta <- apply(tmp,1,mean,na.rm=T)


	# run nonIterative model
	lambda <- calcRangeLambda_dyn(delta,obs,stepsize=stepsize,deltaPk=deltaPk,deltaPt=deltaPt,deltaPkt=deltaPkt)
	
  print("all lambdas")
  print(lambda)


	edges_all_lpSolve = vector()
	baseline_all_lpSolve = vector()
	edges_all_gurobi = vector()
	baseline_all_gurobi = vector()
	sumTime_lpSolve = 0
	sumTime_gurobi = 0
	
	for(lamd in lambda)
	{
	
		print (paste("current lambda", lamd))
		time_lpSolve = system.time((res_lpSolve = LPfunction_lpSolve(obs=obs, delta=delta,lambda=lamd,b=bvec,n=n,K=K,T_=T_,annot,prior=NULL,sourceNode=NULL,sinkNode=NULL,all.int=FALSE,all.pos=FALSE,deltaPk=deltaPk,deltaPt=deltaPt,deltaPkt=deltaPkt)))

		sumTime_lpSolve = sumTime_lpSolve + time_lpSolve
		adja <- getAdja(res_lpSolve$res,n)
		baseline <- getBaseline(res_lpSolve$res,n=n)
		edges_all_lpSolve <- rbind(edges_all_lpSolve,c(t(adja)))
		baseline_all_lpSolve <- rbind(baseline_all_lpSolve, baseline)
		
		time_gurobi = system.time((res_gurobi = LPfunction_gurobi(obs=obs, delta=delta,lambda=lamd,b=bvec,n=n,K=K,T_=T_,annot,prior=NULL,sourceNode=NULL,sinkNode=NULL,all.int=FALSE,all.pos=FALSE,deltaPk=deltaPk,deltaPt=deltaPt,deltaPkt=deltaPkt)))
		
		sumTime_gurobi = sumTime_gurobi + time_gurobi
		
		adja <- getAdja_gurobi(res_gurobi$res,n, res_gurobi$annot)
		baseline <- getBaseline_gurobi(res_gurobi$res,n=n)
		edges_all_gurobi <- rbind(edges_all_gurobi,c(t(adja)))
		baseline_all_gurobi <- rbind(baseline_all_gurobi, baseline)



		checkEquals(res_lpSolve$res$solution, res_gurobi$res$x)
		checkEquals(res_lpSolve$res$objval, res_gurobi$res$objval)
#		checkEquals(edges_all_lpSolve, edges_all_gurobi)
#		checkEquals(baseline_all_lpSolve, baseline_all_gurobi)
	}
	
	write(c(sumTime_gurobi, sumTime_lpSolve), file="execution time - lpSolve_gurobi", sep=" ")
	print (sumTime_gurobi)
	print (sumTime_lpSolve)
}
