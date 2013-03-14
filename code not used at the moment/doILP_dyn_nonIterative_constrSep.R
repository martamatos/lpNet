doILP_dyn_nonIterative_constrSep <-
function(obs,delta,lambda,lambda2,b,n,K,T_,annot,previousNet=NULL,baseline=NULL,previousBaseline=NULL,prior=NULL,sourceNode=NULL,sinkNode=NULL,all.int=FALSE,all.pos=FALSE)
{

	# for now K is the numebr of time points
	annot <- getEdgeAnnot(n=n,T_=T_)
  if(all.int){
		delta <- rep(1,n)
  }

	nConstr = n*K*T_ + 2*n*n*K*(T_-1)
	
  ## weight matrix 
  if(all.pos) W <- matrix(0,nrow=nConstr,ncol=T_*n*n+n)
  else W <- matrix(0,nrow=nConstr,ncol=2*T_*n*n+n)
  colnames(W) <- annot
  # direction of inequation
  f.dir <- rep("<=",nConstr)

  ## convert observations into matrix-format
  # Vector of numeric values for the right-hand sides of the constraints
  bvec <- rep(0,nConstr)
  J <- seq(1,n)
  count <- 1
  slack_var <- rep(FALSE,nConstr) # TRUE if the sign of equation is changed
  
  
  # include code for t=1
  
   for (t in 1:1)
	 {
		for(k in 1:K)
		{
			for(i in 1:n)
			{
				# if the entry in b is 1 then the gene is active
				if(b[(k-1)*n + i]==1){
				# if the observation=NA, just do nothing
					if(!is.na(obs[i,k,t])){
#						print(paste("first t ", i, sep=""))
						# if observation of gene i after knockdown k is acitve
						if(obs[i,k,t]>= delta[i])
						{
							slack_var[count] <- TRUE # what's this variable for?
							# set offset parameter (baseline of gene i)
							if(all.pos)
							{
								W[count,i+(T_*n*n)] <- 1 
								# sum
								for(j in J[J!=i])
								{
									# positive parameter
									id <- which(annot==paste("w+",j,i,t,sep="_"))
									W[count,id] <- obs[j,k,t]
								}
							}
							else
							{
								W[count,i+(2*T_*n*n)] <- 1  
								for(j in J[J!=i])
								{
									# positive parameter
									id <- which(annot==paste("w+",j,i,t,sep="_"))
									W[count,id] <- obs[j,k,t]
									# negative parameter
									id <- which(annot==paste("w-",j,i,t,sep="_"))
									W[count,id] <- -obs[j,k,t]
								} 
							}
							f.dir[count] <- ">="
							bvec[count] <- delta[i]
						}
						# if observation of gene i after knockdown k is NOT acitve
						if(obs[i,k,t]< delta[i])
						{
							# set offset parameter (baseline of gene i)
							if(all.pos)
							{
								W[count,i+(T_*n*n)] <- 1		
								# sum
								for(j in J[J!=i])
								{
									# positive parameter
									id <- which(annot==paste("w+",j,i,t,sep="_"))
									W[count,id] <- obs[j,k,t]
								}
							}
							else
							{
								W[count,i+(2*T_*n*n)] <- 1		
								# sum
								for(j in J[J!=i])
								{
									# positive parameter
									id <- which(annot==paste("w+",j,i,t,sep="_"))
									W[count,id] <- obs[j,k,t]
									# negative parameter
									id <- which(annot==paste("w-",j,i,t,sep="_"))
									W[count,id] <- -obs[j,k,t]
								}
							}
							f.dir[count] <- "<="
							bvec[count] <- 0
						}
					}
				}
				count <- count+1	
			}	# end i		
		} # end k
	} # end t

  
  
  for (t in 2:T_)
  {
		for(k in 1:K)
		{
			for(i in 1:n)
			{
				# if the entry in b is 1 then the gene is active
				if(b[(k-1)*n + i]==1){
				# if the observation=NA, just do nothing
					if(!is.na(obs[i,k,t])){
						# if observation of gene i after knockdown k is acitve
						if(obs[i,k,t]>= delta[i])
						{
#						print("in")
#						print(t)
#						print(i)
#						print(count)
							slack_var[count] <- TRUE # what's this variable for?
							# set offset parameter (baseline of gene i)
							if(all.pos)
							{
								W[count,i+(T_*n*n)] <- 1 
								# sum
								for(j in J[J!=i])
								{
									# positive parameter
									id <- which(annot==paste("w+",j,i,t,sep="_"))
									W[count,id] <- obs[j,k,t]
								}
							}
							else
							{
								W[count,i+(2*T_*n*n)] <- 1  
								for(j in J[J!=i])
								{
									# positive parameter
									id <- which(annot==paste("w+",j,i,t,sep="_"))
									W[count,id] <- obs[j,k,t]
									# negative parameter
									id <- which(annot==paste("w-",j,i,t,sep="_"))
									W[count,id] <- -obs[j,k,t]
								} 
							}
							f.dir[count] <- ">="
							bvec[count] <- delta[i]
						}
						# if observation of gene i after knockdown k is NOT acitve
						if(obs[i,k,t]< delta[i])
						{
							# set offset parameter (baseline of gene i)
							if(all.pos)
							{
								W[count,i+(T_*n*n)] <- 1		
								# sum
								for(j in J[J!=i])
								{
									# positive parameter
									id <- which(annot==paste("w+",j,i,t,sep="_"))
									W[count,id] <- obs[j,k,t]
								}
							}
							else
							{
								W[count,i+(2*T_*n*n)] <- 1		
								# sum
								for(j in J[J!=i])
								{
									# positive parameter
									id <- which(annot==paste("w+",j,i,t,sep="_"))
									W[count,id] <- obs[j,k,t]
									# negative parameter
									id <- which(annot==paste("w-",j,i,t,sep="_"))
									W[count,id] <- -obs[j,k,t]
								}
							}
							f.dir[count] <- "<="
							bvec[count] <- 0
						}
					} # end if
				} # end if
				count = count + 1
			} # end i
						
						
						
			# positive edges 
			for(i in (n+1):(n*n+n))
			{
				ii <- (i - 1) %/% n
				jj <- (i - 1) %% n + 1
				
				# if the entry in b is 1 then the gene is active
				if(b[(k-1)*n + ii]==1){
					# if the observation=NA, just do nothing
					if(!is.na(obs[ii,k,t])){
						if(all.pos){
							if (ii != jj)
							{
								# positive parameter
								id <- which(annot==paste("w+",jj,ii,t-1,sep="_"))
								W[count,id] <- 1 #obs[j,k,t]
								id <- which(annot==paste("w+",jj,ii,t,sep="_"))
								W[count,id] <- 1 #obs[j,k,t]
							}
						}
						else{
							if (ii != jj)
							{
								# positive parameter
								id <- which(annot==paste("w+",jj,ii,t-1,sep="_"))
								W[count,id] <- 1 #obs[j,k,t]
								id <- which(annot==paste("w+",jj,ii,t,sep="_"))
								W[count,id] <- 1 #obs[j,k,t]
							}
						}
						f.dir[count] <- "="
						bvec[count] <- 0
					}
				}
				count = count + 1
			} # end i
	
			
			# negative edges
			for(i in (n+1):(n*n+n))
			{
				ii <- (i - 1) %/% n
				jj <- (i - 1) %% n + 1
				
				# if the entry in b is 1 then the gene is active
				if(b[(k-1)*n + ii]==1){
					# if the observation=NA, just do nothing
					if(!is.na(obs[ii,k,t])){
						if(!all.pos){
							if (ii != jj)
							{
								# negative parameter
								id <- which(annot==paste("w-",jj,ii,t-1,sep="_"))
								W[count,id] <- -1 #obs[j,k,t]
								id <- which(annot==paste("w-",jj,ii,t,sep="_"))
								W[count,id] <- -1 #obs[j,k,t]
							}
						}
						f.dir[count] <- "="
						bvec[count] <- 0
					}
				}
				count = count + 1
			} # end i
			
		}	# end k
	} # end t

	
#  print(W)



#	print(f.dir)
	if(lambda!=0){
		# add slack variables
		sl <- matrix(0,nrow=nConstr,ncol=nConstr)
		annot_s <- paste("s",seq(1,nConstr*2),sep="_")
		colnames(sl) <- annot_s[1:nConstr]

		# attention: for the constr. where observation is smaller than threshold
		xi <- vector()
		for(j in 1:length(f.dir)){
			if(f.dir[j]==">=") xi <- c(xi,0)
			if(f.dir[j]=="<=") xi <- c(xi,-1)
			if(f.dir[j]=="=") xi = c(xi,1)
		}
		diag(sl) <- xi
		W <- cbind(W,sl)
		
		sl <- matrix(0,nrow=nConstr,ncol=nConstr)
		colnames(sl) <- annot_s[(nConstr+1):(2*nConstr)]
	
		xi = c()
		for(j in 1:length(f.dir)){
			if(f.dir[j]==">=") xi <- c(xi,0)
			if(f.dir[j]=="<=") xi <- c(xi,0)
			if(f.dir[j]=="=") xi = c(xi,-1)
		}

		diag(sl) <- xi
		W <- cbind(W,sl)


		#build objective function
		if(all.pos){
			cvec <- c(rep(1,T_*n*n),rep(1,n),rep(1/lambda,length(annot_s)))
			# self-activation is not allowed
			id_self <- c(which(annot==paste("w+",seq(1,n),seq(1,n),seq(1,T_),sep="_")))
		}
		else{
			cvec <- c(rep(1,T_*n*n),rep(1,T_*n*n),rep(1,n),rep(1/lambda,length(annot_s)))
			# self-activation is not allowed
			id_self <- c(which(annot==paste("w+",seq(1,n),seq(1,n),seq(1,T_),sep="_")),which(annot==paste("w-",seq(1,n),seq(1,n),seq(1,T_),sep="_")))
		}
		cvec[id_self] <- 0
		names(cvec) <- c(annot,annot_s)
  }
  else {	#build objective function with no slack variable
		if(all.pos){
			cvec <- c(rep(1,T_*n*n),rep(1,n)) 
			# self-activation is not allowed
			id_self <- c(which(annot==paste("w+",seq(1,n),seq(1,n),seq(1,T_),sep="_")))
		}
		else{
			cvec <- c(rep(1,T_*n*n),rep(1,T_*n*n),rep(1,n)) 
			# self-activation is not allowed
			id_self <- c(which(annot==paste("w+",seq(1,n),seq(1,n),seq(1,T_),sep="_")),which(annot==paste("w-",seq(1,n),seq(1,n),seq(1,T_),sep="_")))
		}
		cvec[id_self] <- 0
		names(cvec) <- c(annot)
  }
 
#	print("W on LP_dyn")
#	print(W)
#	print("bvec")
#	print(bvec)
#	print("cvec")
#	print(cvec)

  ## if there is a prior
  if(!is.null(prior))
  {
		for(i in 1:length(prior))
		{
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



