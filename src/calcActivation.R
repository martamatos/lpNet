calcActivation <-
function(T_nw,b,n,K)
{
  kds <- matrix(b,nrow=n,ncol=K)
  activation_mat <- matrix(NA,nrow=n,ncol=K)
  
  for(k in 1:dim(kds)[2]){
		nw <- T_nw
		inflow <- vector("list", length=n)
		# in-degree for each edge:
		in_deg <- apply(abs(nw),2,sum)
		# find root_nodes
		root_nodes <- which(in_deg == 0)
		# if no root nodes: stop
		if(length(root_nodes)==0){
			cat("Error: there are no root nodes\n")
		}
		done <- vector()
		# process root_nodes: then delete there outgoing edges: set new indegree
		for(i in root_nodes){
			# root nodes are inactive if kd 
			if(kds[i,k]==0){
				inflow[[i]] <- c(inflow[[i]],0)
			}
			# otherwise they are active and its children (if not kd)
			else{
				inflow[[i]] <- c(inflow[[i]],1)
				children <- which(nw[i,]!=0)
				for(c in children){
					# if children kd -> they are inactive
					if(kds[c,k]==0){
						inflow[[c]] <- c(inflow[[c]],0)
						nw[c,] <- 0
						done <- c(done,c)
					}
					else{
						if(nw[i,c]>0) inflow[[c]] <- c(inflow[[c]],1)
						else inflow[[c]] <- c(inflow[[c]],-1)
						in_deg[c] <- in_deg[c]-1
					}
				}
			}
			nw[i,] <- 0
			in_deg <- apply(abs(nw),2,sum)
			done <- c(done,i)
		}
		# now proceed with nodes where in_deg is zero and which are not done already
		ids_tmp <- which(in_deg==0)
		ids <- ids_tmp[!ids_tmp %in% done]
		
		while(length(ids)>0){
		
			for(i in ids){
				parents <- which(T_nw[,i]!=0)
				
				for(pa in parents){
					if(sum(inflow[[pa]])<=0) inflow[[i]] <- c(inflow[[i]],0)
					else{
						if(T_nw[pa,i]>0) inflow[[i]] <- c(inflow[[i]],1)
						else inflow[[i]] <- c(inflow[[i]],-1)
					}
				}
				children <- which(nw[i,]!=0)
				
				for(c in children){
					# if children kd -> they are inactive
					if(kds[c,k]==0){
						inflow[[c]] <- c(inflow[[c]],0)
						nw[c,] <- 0
						in_deg <- apply(abs(nw),2,sum)
						done <- c(done,c)
					}
					else{
						if(sum(inflow[[i]])<=0) inflow[[c]] <- c(inflow[[c]],0)
						else{
							if(nw[i,c]>0) inflow[[c]] <- c(inflow[[c]],1)
							else inflow[[c]] <- c(inflow[[c]],-1)
						}
						in_deg[c] <- in_deg[c]-1
					}
				}
				done <- c(done,i)
				ids_tmp <- which(in_deg==0)
				ids <- ids_tmp[!ids_tmp %in% done]
			}
		}
		# if no indegree is zero and some nodes are still undone: there is a loop
		undone <- which(!seq(1,n)%in%done)
		incom <- unlist(lapply(inflow,sum))
		
		while(length(undone)>0){
			# if any undone node which is not root has inflow>0 start there
			ids <- undone[incom[undone]!=0]
			if(length(ids)>0){
				for(i in ids){
					# if node is kd: inflow is zero (otherwise just let inflow like it is
					if(kds[i,k]==0){
						inflow[[i]] <- c(inflow[[i]],0)
					}
					else{
						children <- which(nw[i,]!=0)
						for(c in children){
							# if children kd -> they are inactive
							if(kds[c,k]==0){
							inflow[[c]] <- c(inflow[[c]],0)
							}
							else{
								if(sum(inflow[[i]])<=0) inflow[[c]] <- c(inflow[[c]],0)
								else{
									if(nw[i,c]>0) inflow[[c]] <- c(inflow[[c]],1)
									else inflow[[c]] <- c(inflow[[c]],-1)
								}
							}
						}
					}
					done <- c(done,i)
				}
			}
			# if no node is active: rest is inactive
			else{
				for(i in undone){
					inflow[[i]] <- c(inflow[[i]],0)
					done <- c(done,i)
				}
			}
			undone <- which(!seq(1,n)%in%done)
			incom <- unlist(lapply(inflow,sum))
		}
		tmp <- unlist(lapply(inflow,sum))
		activation_mat[,k] <- apply(cbind(rep(0,n),tmp),1,max)
	}
  activation_mat[activation_mat!=0] <- 1

  return(activation_mat)
}




#
# calculates only the effects of knockdowns, assuming only pos edges. 
#  (so that there are no inactive nodes due to inhibiting edges)
#

calcActivation_dyn <- function(T_nw,b,n,K)
{

  kds <- matrix(b,nrow=n,ncol=K)
  activation_mat <- matrix(NA,nrow=n,ncol=K)
  
  for(k in 1:dim(kds)[2]){
		nw <- T_nw
		inflow <- vector("list", length=n)
		# in-degree for each edge:
		in_deg <- apply(abs(nw),2,sum)
		# find root_nodes
		root_nodes <- which(in_deg == 0)
		# if no root nodes: stop
		if(length(root_nodes)==0){
			cat("Error: there are no root nodes\n")
		}
		done <- vector()
		# process root_nodes: then delete there outgoing edges: set new indegree
		for(i in root_nodes){
			# root nodes are inactive if kd 
			if(kds[i,k]==0){
				inflow[[i]] <- c(inflow[[i]],0)
			}
			# otherwise they are active and its children (if not kd)
			else{
				inflow[[i]] <- c(inflow[[i]],1)
				children <- which(nw[i,]!=0)
				for(c in children){
					# if children kd -> they are inactive
					if(kds[c,k]==0){
						inflow[[c]] <- c(inflow[[c]],0)
						nw[c,] <- 0
						done <- c(done,c)
					}
					else{
						if(nw[i,c]!=0) inflow[[c]] <- c(inflow[[c]],1)
						in_deg[c] <- in_deg[c]-1
					}
				}
			}
			nw[i,] <- 0
			in_deg <- apply(abs(nw),2,sum)
			done <- c(done,i)
		}
		# now proceed with nodes where in_deg is zero and which are not done already
		ids_tmp <- which(in_deg==0)
		ids <- ids_tmp[!ids_tmp %in% done]
		
		while(length(ids)>0){
		
			for(i in ids){
				parents <- which(T_nw[,i]!=0)
				
				for(pa in parents){
					if(sum(inflow[[pa]])<=0) inflow[[i]] <- c(inflow[[i]],0)
					else{
						if(T_nw[pa,i]!=0) inflow[[i]] <- c(inflow[[i]],1)
					}
				}
				children <- which(nw[i,]!=0)
				
				for(c in children){
					# if children kd -> they are inactive
					if(kds[c,k]==0){
						inflow[[c]] <- c(inflow[[c]],0)
						nw[c,] <- 0
						in_deg <- apply(abs(nw),2,sum)
						done <- c(done,c)
					}
					else{
						if(sum(inflow[[i]])<=0) inflow[[c]] <- c(inflow[[c]],0)
						else{
							if(nw[i,c]!=0) inflow[[c]] <- c(inflow[[c]],1)
						}
						in_deg[c] <- in_deg[c]-1
					}
				}
				done <- c(done,i)
				ids_tmp <- which(in_deg==0)
				ids <- ids_tmp[!ids_tmp %in% done]
			}
		}
		# if no indegree is zero and some nodes are still undone: there is a loop
		undone <- which(!seq(1,n)%in%done)
		incom <- unlist(lapply(inflow,sum))
		
		while(length(undone)>0){
			# if any undone node which is not root has inflow>0 start there
			ids <- undone[incom[undone]!=0]
			if(length(ids)>0){
				for(i in ids){
					# if node is kd: inflow is zero (otherwise just let inflow like it is
					if(kds[i,k]==0){
						inflow[[i]] <- c(inflow[[i]],0)
					}
					else{
						children <- which(nw[i,]!=0)
						for(c in children){
							# if children kd -> they are inactive
							if(kds[c,k]==0){
								inflow[[c]] <- c(inflow[[c]],0)
							}
							else{
								if(sum(inflow[[i]])<=0) inflow[[c]] <- c(inflow[[c]],0)
								else{
									if(nw[i,c]!=0) inflow[[c]] <- c(inflow[[c]],1)
								}
							}
						}
					}
					done <- c(done,i)
				}
			}
			# if no node is active: rest is inactive
			else{
				for(i in undone){
					inflow[[i]] <- c(inflow[[i]],0)
					done <- c(done,i)
				}
			}
			undone <- which(!seq(1,n)%in%done)
			incom <- unlist(lapply(inflow,sum))

		}
		tmp <- unlist(lapply(inflow,sum))
		activation_mat[,k] <- apply(cbind(rep(0,n),tmp),1,max)

	}
  activation_mat[activation_mat!=0] <- 1

  return(activation_mat)
}
