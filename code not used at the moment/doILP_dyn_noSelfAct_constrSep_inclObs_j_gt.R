doILP_dyn_noSelfAct_constrSep_inclObs_j_gt <-
function(obs,delta,lambda,lambda2,b,n,K,T_=NULL,annot,previousNet=NULL,baseline=NULL,previousBaseline=NULL,prior=NULL,sourceNode=NULL,sinkNode=NULL,all.int=FALSE,all.pos=FALSE)
{
	print(paste("lambda2 beg: ", lambda2, sep=""))
	# for now K is the numebr of time points

  if(all.int){
		delta <- rep(1,n)
  }

	nConstr = K*(n+n*n)
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
		}
		#print(count)
		
		# add "dynamic" constraints
		for(i in (n+1):(n*n+n)){

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
								W[count,id] <- obs[jj,k]
							}
						}
						else{
							if (ii != jj){
								# positive parameter
								id <- which(annot==paste("w+",jj,ii,sep="_"))
								W[count,id] <- obs[jj,k]
								
								# negative parameter
								id <- which(annot==paste("w-",jj,ii,sep="_"))
								W[count,id] <- obs[jj,k]
							} 
						}
					f.dir[count] <- ">="
					bvec[count] <- abs(previousNet[jj,ii])
				} # end if
			} # end if
			count <- count+1
		} # end i
  } # end k
  
  
  
#  print("fdir")
#  print(f.dir)
  #print(W)
#  lambda2=lambda
  ## now add slack varibles to W: 
  
	sl <- matrix(0,nrow=nConstr,ncol=nConstr)
	count_slack = 1
	annot_s = vector()
	
	# attention: for the constr. where observation is smaller than threshold
	xi <- vector()
#		fdir_len = length(f.dir)
	for(j in 1:nConstr){
		if(f.dir[j]==">=" && j <= n){
			xi <- c(xi,0)
			annot_s = c(annot_s, paste("stI", count_slack, sep="_"))
		}
		if(f.dir[j]=="<=" && j <= n){
			xi <- c(xi,-1)
			annot_s = c(annot_s, paste("stI", count_slack, sep="_"))
		}
		if(f.dir[j]==">=" && j > n){
			xi = c(xi,1)
			annot_s = c(annot_s, paste("dyn", count_slack, sep="_"))
		}
		if(f.dir[j]=="<=" && j > n){
			xi = c(xi,222)
			annot_s = c(annot_s, paste("dyn", count_slack, sep="_"))
		}
		count_slack = count_slack + 1
	}

	colnames(sl) <- annot_s[1:nConstr]
	diag(sl) <- xi
	W <- cbind(W,sl)
	
	sl <- matrix(0,nrow=nConstr,ncol=nConstr)

	xi = c()
	for(j in 1:nConstr){
		if(j <= n){
			xi <- c(xi,111)
			annot_s = c(annot_s, paste("stNA", count_slack, sep="_"))
		}
		if(f.dir[j]==">=" && j > n){
			xi = c(xi,-1)
			annot_s = c(annot_s, paste("dyn", count_slack, sep="_"))
		}
		if(f.dir[j]=="<=" && j > n){
			xi = c(xi,222)
			annot_s = c(annot_s, paste("dyn", count_slack, sep="_"))
		}
		count_slack = count_slack + 1
	}

	colnames(sl) <- annot_s[(nConstr+1):(2*nConstr)]
	diag(sl) <- xi
	W <- cbind(W,sl)
	
	
	
	cols_dyn = grep("dyn",colnames(W), fixed=TRUE)
	cols_dyn = matrix(c(cols_dyn), nrow=1, ncol=length(cols_dyn))
	colnames(cols_dyn) = colnames(W[,cols_dyn])
	
	rows_dyn = which(W[,cols_dyn]!=0,arr.ind=T )[,1]
	entries_dyn = matrix(c(rows_dyn, cols_dyn), nrow=length(rows_dyn), ncol=2, byrow=F)
	
	W[which(W[,cols_dyn]==222,arr.ind=T)[,1], cols_dyn]=0

	
	cols_staticI = grep("stI", colnames(W), fixed=TRUE)
	cols_staticI = matrix(c(cols_staticI), nrow=1, ncol=length(cols_staticI))
	colnames(cols_staticI) = colnames(W[,cols_staticI])

	
	cols_staticNA = grep("stNA", colnames(W), fixed=TRUE)
	cols_staticNA = matrix(c(cols_staticNA), nrow=1, ncol=length(cols_staticNA))
	colnames(cols_staticNA) = colnames(W[,cols_staticNA]) 
	
	rows_static = which(W[,cols_staticNA]==111,arr.ind=T)[,1]

	entries_staticI = matrix(c(rows_static, cols_staticI), nrow=length(rows_static), ncol=2, byrow=F)
	entries_staticNA = matrix(c(rows_static, cols_staticNA), nrow=length(rows_static), ncol=2, byrow=F)

	W[entries_staticNA] = 0
	
		
	if(lambda!=0 && lambda2!=0)
	{
		#build objective function with lambda and lambda2
		if(all.pos){
			cvec <- c(rep(1,n*n),rep(1,n),
								c(1/lambda*abs(W[entries_staticI])),
								c(1/lambda*abs(W[entries_staticNA])),
								c(1/lambda2*abs(W[entries_dyn])))
			# self-activation is not allowed
			id_self <- c(which(annot==paste("w+",seq(1,n),seq(1,n),seq(1,T_),sep="_")))
		}
		else{
			cvec <- c(rep(1,n*n),rep(1,n*n),rep(1,n),
								c(1/lambda*abs(W[entries_staticI])),
								c(1/lambda*abs(W[entries_staticNA])),
								c(1/lambda2*abs(W[entries_dyn])))
			# self-activation is not allowed
			id_self <- c(which(annot==paste("w+",seq(1,n),seq(1,n),seq(1,T_),sep="_")),which(annot==paste("w-",seq(1,n),seq(1,n),seq(1,T_),sep="_")))
		}
		cvec[id_self] <- 0
		names(cvec) <- c(annot,	colnames(cols_staticI),	colnames(cols_staticNA), colnames(cols_dyn))
  }
	else if(lambda!=0 && lambda2==0)
	{
		#build objective function with lambda
		if(all.pos){
			cvec <- c(rep(1,n*n),rep(1,n),
								c(1/lambda*abs(W[entries_staticI])),
								c(1/lambda*abs(W[entries_staticNA])))
			# self-activation is not allowed
			id_self <- c(which(annot==paste("w+",seq(1,n),seq(1,n),seq(1,T_),sep="_")))
		}
		else{
			cvec <- c(rep(1,n*n),rep(1,n*n),rep(1,n),
								c(1/lambda*abs(W[entries_staticI])),
								c(1/lambda*abs(W[entries_staticNA])))
			# self-activation is not allowed
			id_self <- c(which(annot==paste("w+",seq(1,n),seq(1,n),seq(1,T_),sep="_")),which(annot==paste("w-",seq(1,n),seq(1,n),seq(1,T_),sep="_")))
		}
		cvec[id_self] <- 0
		names(cvec) <- c(annot,	colnames(cols_staticI),	colnames(cols_staticNA))
  }
  else if(lambda==0 && lambda2!=0)
  {
		#build objective function with lambda2
		if(all.pos){
			cvec <- c(rep(1,n*n),rep(1,n),
								c(1/lambda2*abs(W[entries_dyn])))
			# self-activation is not allowed
			id_self <- c(which(annot==paste("w+",seq(1,n),seq(1,n),seq(1,T_),sep="_")))
		}
		else{
			cvec <- c(rep(1,n*n),rep(1,n*n),rep(1,n),
								c(1/lambda2*abs(W[entries_dyn])))
			# self-activation is not allowed
			id_self <- c(which(annot==paste("w+",seq(1,n),seq(1,n),seq(1,T_),sep="_")),which(annot==paste("w-",seq(1,n),seq(1,n),seq(1,T_),sep="_")))
		}
		cvec[id_self] <- 0
		names(cvec) <- c(annot, colnames(cols_dyn))
  }
  else if(lambda==0 && lambda2==0) 
  {	#build objective functio with no slack variables
		if(all.pos){
			cvec <- c(rep(1,n*n),rep(1,n)) 
			# self-activation is not allowed
			id_self <- c(which(annot==paste("w+",seq(1,n),seq(1,n),seq(1,T_),sep="_")))
		}
		else{
			cvec <- c(rep(1,n*n),rep(1,n*n),rep(1,n)) 
			# self-activation is not allowed
			id_self <- c(which(annot==paste("w+",seq(1,n),seq(1,n),seq(1,T_),sep="_")),which(annot==paste("w-",seq(1,n),seq(1,n),seq(1,T_),sep="_")))
		}
		cvec[id_self] <- 0
		names(cvec) <- c(annot)
  }
 
 
#	print("W on LP_dyn")
#	print(W)
#	print("lambda2")
#	print(lambda2)
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



