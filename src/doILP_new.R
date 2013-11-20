doILP_new <-
function(obs,delta,lambda,b,n,K,T_=NULL,annot,prior=NULL,sourceNode=NULL,sinkNode=NULL,all.int=FALSE,all.pos=FALSE,deltaPk=FALSE,deltaPt=FALSE,deltaPkt=FALSE)
{
  if(all.int){
		delta <- rep(1,n)
  }

  ## weight matrix of dim ((K*n)x(2n²+n)) (w_i^0)
  if(all.pos) W <- matrix(0,nrow=K*n,ncol=n*n+n)
  else W <- matrix(0,nrow=K*n,ncol=2*n*n+n)
  colnames(W) <- annot
  # direction of inequation
  f.dir <- rep("<=",n*K)
  
  ## convert observations into matrix-format
  # Vector of numeric values for the right-hand sides of the constraints
  bvec <- rep(0,n*K)
  J <- seq(1,n)
  count <- 1
  slack_var <- rep(FALSE,K*n) # TRUE if the sign of equation is changed
  
  if (deltaPk==F){
		if(all.int){
			delta <- rep(1,n)
		}
		
		for(k in 1:K){
			for(i in 1:n){
				# if the entry in b is 1 then gene count has been silenced and thus is not activated after knockdown k
				if(b[count]==1){
				# if the observation=NA, just do nothing
				if(!is.na(obs[i,k])){
					# if observation of gene i after knockdown k is acitve
					if(obs[i,k]>= delta[i]){
						slack_var[count] <- TRUE
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
		}
	}
  else if (deltaPk==T){
		if(all.int){
			delta <- matrix(rep(1,n*K), nrow=n, ncol=K)
		}
		
		for(k in 1:K){
			for(i in 1:n){
				# if the entry in b is 1 then gene count has been silenced and thus is not activated after knockdown k
				if(b[count]==1){
				# if the observation=NA, just do nothing
				if(!is.na(obs[i,k])){
					# if observation of gene i after knockdown k is acitve
					if(obs[i,k]>= delta[i,k]){
						slack_var[count] <- TRUE
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
						bvec[count] <- delta[i,k]
					}
					# if observation of gene i after knockdown k is NOT acitve
					if(obs[i,k]< delta[i,k]){
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
		}
	}
  
  ## now add slack varibles to W: 
  #slack variables=length(b)
  if(lambda!=0){
		sl <- matrix(0,nrow=dim(W)[1],ncol=dim(W)[1])
		annot_s <- paste("s",seq(1,dim(W)[1]),sep="_")
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
			cvec <- c(rep(1,n*n),rep(1,n),rep(1/lambda,dim(sl)[1]))
			# self-activation is not allowed
			id_self <- c(which(annot==paste("w+",seq(1,n),seq(1,n),sep="_")))
		}
		else{
			cvec <- c(rep(1,n*n),rep(1,n*n),rep(1,n),rep(1/lambda,dim(sl)[1]))
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
  

  ## condition that each node which is not End hast at least delta[i] outgoing edges
  if(!is.null(sinkNode)){
	W_tmp1 <- vector()
	gene_tmp <- seq(1,n)[-sinkNode]
	for(i in gene_tmp){
	  # outgoing edge can come from all nodes except itself
	  tmp <- seq(1,n)[-i]
	  if(length(tmp)>1){
			# for negative and positive parameter
			annot_pos <- paste("w+",i,tmp,sep="_")
			
			if(!all.pos) annot_neg <- paste("w-",i,tmp,sep="_")
			add_row <- rep(0,length(cvec))
			add_row[which(annot%in%annot_pos)] <- 1
			
			if(!all.pos) add_row[which(annot%in%annot_neg)] <- 1
			W_tmp1 <- rbind(W_tmp1,as.double(add_row))
			bvec <- c(bvec,delta[i])
			f.dir <- c(f.dir,">=")
	  }
	}
	W <- rbind(W,W_tmp1)
  }
  
  ## conditions that each node which is not Start has at least delta[i] incoming edges
  if(!is.null(sourceNode)){
		W_tmp2 <- vector()
		gene_tmp <- seq(1,n)[-sourceNode]
		for(i in gene_tmp){
			# incoming edge can come from all nodes except itself
			tmp <- seq(1,n)[-i]
			if(length(tmp)>1){
				annot_pos <- paste("w+",tmp,i,sep="_")
				
				if(!all.pos) annot_neg <- paste("w-",tmp,i,sep="_")
				add_row <- rep(0,length(cvec))
				add_row[which(annot%in%annot_pos)] <- 1
				
				if(!all.pos) add_row[which(annot%in%annot_neg)] <- 1
				W_tmp2 <- rbind(W_tmp2,as.double(add_row))
				bvec <- c(bvec,delta[i])
				f.dir <- c(f.dir,">=")
			}
		}
		W <- rbind(W,W_tmp2)
  }

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
  return(res)
}
