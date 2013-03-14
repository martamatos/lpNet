doILP_dyn <-
function(obs,delta,lambda,lambda2,b,n,K,T_,annot,previousNet=NULL,baseline=NULL,previousBaseline=NULL,prior=NULL,sourceNode=NULL,sinkNode=NULL,all.int=FALSE,all.pos=FALSE)
{

	# for now K is the numebr of time points

  if(all.int){
		delta <- rep(1,n)
  }

	nConstr = n*K*(T_-1)
	
  ## weight matrix of dim ((K*n)x(2n²+n)) (w_i^0)
  if(all.pos) W <- matrix(0,nrow=nConstr,ncol=n*n+n)
  else W <- matrix(0,nrow=nConstr,ncol=2*n*n+n)
  colnames(W) <- annot
  # direction of inequation
  f.dir <- rep("<=",nConstr)

  ## convert observations into matrix-format
  # Vector of numeric values for the right-hand sides of the constraints
  bvec <- rep(0,nConstr)
  J <- seq(1,n)
  count <- 1
  slack_var <- rep(FALSE,nConstr) # TRUE if the sign of equation is changed
  
  
  for (t in 2:T_){
		for(k in 1:K){
			for(i in 1:n){
				# if the entry in b is 1 then the gene is active
				if(b[(k-1)*n + i]==1){
					# if the observation=NA, just do nothing
					if(!is.na(obs[i,k,t])){
						# if observation of gene i after knockdown k is acitve
						if(obs[i,k,t]>= delta[i]){
							slack_var[count] <- TRUE # what's this variable for?
							# set offset parameter (baseline of gene i)
							if(all.pos){
								W[count,i+(n*n)] <- 1
								# sum
								for(j in J[J!=i]){
									# positive parameter
									id <- which(annot==paste("w+",j,i,sep="_"))
									W[count,id] <- obs[j,k,t-1]
								}
							}
							else{
								W[count,i+(2*n*n)] <- 1
								for(j in J[J!=i]){
									# positive parameter
									id <- which(annot==paste("w+",j,i,sep="_"))
									W[count,id] <- obs[j,k,t-1]
									# negative parameter
									id <- which(annot==paste("w-",j,i,sep="_"))
									W[count,id] <- -obs[j,k,t-1]
								} 
							}
							f.dir[count] <- ">="
							bvec[count] <- delta[i]
						}
						# if observation of gene i after knockdown k is NOT acitve
						if(obs[i,k,t]< delta[i]){
							# set offset parameter (baseline of gene i)
							if(all.pos){
								W[count,i+(n*n)] <- 1
								# sum
								for(j in J[J!=i]){
									# positive parameter
									id <- which(annot==paste("w+",j,i,sep="_"))
									W[count,id] <- obs[j,k,t-1]
								}
							}
							else{
								W[count,i+(2*n*n)] <- 1
								# sum
								for(j in J[J!=i]){
									# positive parameter
									id <- which(annot==paste("w+",j,i,sep="_"))
									W[count,id] <- obs[j,k,t-1]
									# negative parameter
									id <- which(annot==paste("w-",j,i,sep="_"))
									W[count,id] <- -obs[j,k,t-1]
								}
						}
						f.dir[count] <- "<="
						bvec[count] <- 0
						}
					}
				}
				count <- count+1
			} # end i
		} # end k
	} # end t
		
  
  
#  print("fdir")
#  print(f.dir)
  #print(W)
#  lambda2=lambda
  ## now add slack varibles to W: 
  
	if(lambda!=0){
		sl <- matrix(0,nrow=nConstr,ncol=nConstr)
		annot_s <- paste("s",seq(1,nConstr),sep="_")
		colnames(sl) <- annot_s
		
		# attention: for the constr. where observation is smaller than threshold
		xi <- vector()
		for(j in 1:length(f.dir)){
			if(f.dir[j]==">=") xi <- c(xi,0)
			if(f.dir[j]=="<=") xi <- c(xi,-1)
		}
		diag(sl) <- xi
		W <- cbind(W,sl)
		if(all.pos){
			cvec <- c(rep(1,n*n),rep(1,n),rep(1/lambda,nConstr))
			# self-activation is not allowed
			id_self <- c(which(annot==paste("w+",seq(1,n),seq(1,n),sep="_")))
		}
		else{
			cvec <- c(rep(1,n*n),rep(1,n*n),rep(1,n),rep(1/lambda,nConstr))
			# self-activation is not allowed
			id_self <- c(which(annot==paste("w+",seq(1,n),seq(1,n),sep="_")),which(annot==paste("w-",seq(1,n),seq(1,n),sep="_")))
		}
		cvec[id_self] <- 0
		names(cvec) <- c(annot,annot_s)
	}
	else {
		if(all.pos){
			cvec <- c(rep(1,n*n),rep(1,n)) 
			# self-activation is not allowed
			id_self <- c(which(annot==paste("w+",seq(1,n),seq(1,n),sep="_")))
		}
		else{
			cvec <- c(rep(1,n*n),rep(1,n*n),rep(1,n)) 
			# self-activation is not allowed
			id_self <- c(which(annot==paste("w+",seq(1,n),seq(1,n),sep="_")),which(annot==paste("w-",seq(1,n),seq(1,n),sep="_")))
		}
	cvec[id_self] <- 0
	names(cvec) <- c(annot)
  }
 
 
#	print("W on LP_dyn")
#	print(W)
#	print("cvec")
#	print(cvec)
#	print("bvec")
#	print(bvec)

  ## if there is a prior
  if(!is.null(prior)){
		for(i in 1:length(prior)){
			tmp <- rep(0,dim(W)[2])
			tmp[which(prior[[i]][1]==annot)] <- as.double(prior[[i]][2])
			W <- rbind(W,tmp)
			bvec <- c(bvec,as.double(prior[[i]][4]))
			f.dir <- c(f.dir,prior[[i]][3])
		}
  }
#  print(W)
  ## Maximize the gross margin
  res <- lp("min",cvec,W,f.dir,bvec,all.int=all.int) 
  ## min - direction of optimization
  ## cvec - objective function (Numeric vector of coefficients of objective function)
  ## W - Matrix of numeric constraint coefficients, one row per constraint, one column per variable
  ## f.dir vector of character strings giving the direction of the constraint
  ## bvec - vector of numeric values for the right-hand sides of the constraints


  return(res)
}



