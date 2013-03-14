doILP_dyn_noSelfAct_posNeg <-
function(obs,delta,lambda,b,n,K,T_=NULL,annot,previousNet=NULL,baseline=NULL,previousBaseline=NULL,prior=NULL,sourceNode=NULL,sinkNode=NULL,all.int=FALSE,all.pos=FALSE)
{

	# for now K is the numebr of time points

  if(all.int){
		delta <- rep(1,n)
  }
	nConstr = K*(n+n*n)
  ## weight matrix of dim ((K*n)x(2n²+n)) (w_i^0)
  if(all.pos) W <- matrix(0,nrow=K*(n+n*n),ncol=n*n+n)
  else W <- matrix(0,nrow=K*(n+n*n),ncol=2*n*n+n)
  colnames(W) <- annot
  # direction of inequation
  f.dir <- rep("<=",K*(n+n*n))

  ## convert observations into matrix-format
  # Vector of numeric values for the right-hand sides of the constraints
  bvec <- rep(0,K*(n+n*n))
  J <- seq(1,n)
  count <- 1
  slack_var <- rep(FALSE,K*(n+n*n)) # TRUE if the sign of equation is changed
  
  for(k in 1:K){
		for(i in 1:n){
			# if the entry in b is 1 then the gene is active
			if(b[(k-1)*n + i]==1){
			# if the observation=NA, just do nothing
				
				if(!is.na(obs[i,k])){
					# if observation of gene i after knockdown k is acitve
					if(obs[i,k]>= delta[i]){
						slack_var[count] <- TRUE # what's this variable for?
						# set offset parameter (baseline of gene i)
						if(all.pos){
							W[count,i+(n*n)] <- 1
							# sum
							for(j in J[J!=i]){
								# positive parameter
								id <- which(annot==paste("w+",j,i,sep="_"))
								W[count,id] <- obs[j,k]
							}
						}
						else{
							W[count,i+(2*n*n)] <- 1
							for(j in J[J!=i]){
								# positive parameter
								id <- which(annot==paste("w+",j,i,sep="_"))
								W[count,id] <- obs[j,k]
								# negative parameter
								id <- which(annot==paste("w-",j,i,sep="_"))
								W[count,id] <- -obs[j,k]
								} 
						}
						f.dir[count] <- ">="						
						bvec[count] <- delta[i]					
					}
					# if observation of gene i after knockdown k is NOT acitve
					if(obs[i,k]< delta[i]){
						# set offset parameter (baseline of gene i)
						if(all.pos){
							W[count,i+(n*n)] <- 1
							# sum
							for(j in J[J!=i]){
								# positive parameter
								id <- which(annot==paste("w+",j,i,sep="_"))
								W[count,id] <- obs[j,k]
							}
						}
						else{
							W[count,i+(2*n*n)] <- 1
							# sum
							for(j in J[J!=i]){
								# positive parameter
								id <- which(annot==paste("w+",j,i,sep="_"))
								W[count,id] <- obs[j,k]
								# negative parameter
								id <- which(annot==paste("w-",j,i,sep="_"))
								W[count,id] <- -obs[j,k]
							}
					}
					f.dir[count] <- "<="
					bvec[count] <- 0
					}
				}		
			}
			count <- count+1		
		} # end of i
		#print(count)
		
		for(i in (n+1):(n*n+n))
		{
			ii <- (i - 1) %/% n
			jj <-(i - 1) %% n + 1
			if(b[(k-1)*n + ii]==1){
			# if the observation=NA, just do nothing
				if(!is.na(obs[ii,k])){
					if(all.pos){				
						# sum
						if (ii != jj){
							# positive parameter
							id <- which(annot==paste("w+",jj,ii,sep="_"))
							#"print(count)
							W[count,id] <- 1 #obs[jj,k]
						}
					}
					else{
						if (ii != jj){
							# positive parameter
							id <- which(annot==paste("w+",jj,ii,sep="_"))
							W[count,id] <- 1 #obs[jj,k]
							
							# negative parameter
							id <- which(annot==paste("w-",jj,ii,sep="_"))
							W[count,id] <- -1 #obs[jj,k]
						} 
					}
					f.dir[count] <- ">="
					bvec[count] <- previousNet[jj,ii]
				}
			}
			count <- count+1
		} # end of i
		#print(count)
		
  }
  #print(W)
  ## now add slack varibles to W: 
  #slack variables=length(b)
	if(lambda!=0){
		sl <- matrix(0,nrow=nConstr,ncol=nConstr)
		annot_s <- paste("s",seq(1,nConstr*2),sep="_")
		colnames(sl) <- annot_s[1:nConstr]
#		print(annot_s[1:dim(W)[1]])
		
		# attention: for the constr. where observation is smaller than threshold
		xi <- vector()
		for(j in 1:length(f.dir)){
			if(f.dir[j]==">=" && j <= n) xi <- c(xi,0)
			if(f.dir[j]=="<=" && j <= n) xi <- c(xi,-1)
			if(j > n) xi = c(xi,1)
		}

		diag(sl) <- xi
		W <- cbind(W,sl)
		
		sl <- matrix(0,nrow=nConstr,ncol=nConstr)
		colnames(sl) <- annot_s[(nConstr+1):(2*nConstr)]

		xi = c()

		for(j in 1:length(f.dir)){
			if(f.dir[j]==">=" && j <= n) xi <- c(xi,0)
			if(f.dir[j]=="<=" && j <= n) xi <- c(xi,0)
			if(j > n) xi = c(xi,-1)
		}

		
		diag(sl) <- xi
		W <- cbind(W,sl)
		

		
		
		if(all.pos){
			cvec <- c(rep(1,n*n),rep(1,n),rep(1/lambda,length(annot_s)))
			# self-activation is not allowed
			id_self <- c(which(annot==paste("w+",seq(1,n),seq(1,n),sep="_")))
		}
		else{
			cvec <- c(rep(1,n*n),rep(1,n*n),rep(1,n),rep(1/lambda,length(annot_s)))
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
  
  ## Maximize the gross margin
  res <- lp("min",cvec,W,f.dir,bvec,all.int=all.int) 
  ## min - direction of optimization
  ## cvec - objective function (Numeric vector of coefficients of objective function)
  ## W - Matrix of numeric constraint coefficients, one row per constraint, one column per variable
  ## f.dir vector of character strings giving the direction of the constraint
  ## bvec - vector of numeric values for the right-hand sides of the constraints


  return(res)
}



