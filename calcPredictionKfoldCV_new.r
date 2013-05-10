calcPredictionKfoldCV_new <-function(b,n,K,adja,baseline,obs,delta,rem_entries,active_mu,active_sd,inactive_mu,inactive_sd,muPgene=FALSE,muPgk=FALSE)
{
	# activation matrix is the same regardless of time point
	act_mat <- calcActivation(adja,b,n,K)
	inact_entries = which(act_mat==0, arr.ind=T) # returns (i,k)
	predict = obs
	
	if (muPgene==F & muPgk==F & muPgt==F & muPgkt==F){
	
		for (ent in 1:dim(rem_entries)[1]){
			rem_gene=rem_entries[ent,1]
			rem_k=rem_entries[ent,2]
			
			res=vecMatch(rem_entries[ent,1:2], inact_entries)
			
			# if the removed entry is an inactive node due to some knockdown, then predict as inactivet=rem_entries[ent,3]
			if (any(res)==TRUE){
				predict[rem_gene,rem_k] <- rnorm(1,inactive_mu,inactive_sd)
			}
			else{
				pa <- which(adja[,rem_gene]!=0)
				# if there are no parents: rem_gene is root node and thus is considered to be active since it is not silenced
				if (length(pa) == 0){
					predict[rem_gene,rem_k]  <- rnorm(1,active_mu,active_sd) 
				}
				else{
					in_flow <- baseline[rem_gene]
					
					for(j in 1:length(pa)){
						if (!is.na(obs[pa[j],rem_k]) &  (obs[pa[j],rem_k]>= delta[pa[j]])){
							in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k]),na.rm=T)
						}
					}
					
					if(in_flow >= delta[rem_gene]){
						predict[rem_gene,rem_k]  <- rnorm(1,active_mu,active_sd)
					}
					else{
						predict[rem_gene,rem_k]  <- rnorm(1,inactive_mu,inactive_sd)
					}
				}
			}
		}
	}
	
	else if (muPgene==T){
	
		for (ent in 1:dim(rem_entries)[1]){
			rem_gene=rem_entries[ent,1]
			rem_k=rem_entries[ent,2]
			
			res=vecMatch(rem_entries[ent,1:2], inact_entries)
			
			# if the removed entry is an inactive node due to some knockdown, then predict as inactivet=rem_entries[ent,3]
			if (any(res)==TRUE){
				predict[rem_gene,rem_k] <- rnorm(1,inactive_mu[rem_gene],inactive_sd[rem_gene])
			}
			else{
				pa <- which(adja[,rem_gene]!=0)
				# if there are no parents: rem_gene is root node and thus is considered to be active since it is not silenced
				if (length(pa) == 0){
					predict[rem_gene,rem_k]  <- rnorm(1,active_mu[rem_gene],active_sd[rem_gene]) 
				}
				else{
					in_flow <- baseline[rem_gene]
					
					for(j in 1:length(pa)){
						if (!is.na(obs[pa[j],rem_k]) &  (obs[pa[j],rem_k]>= delta[pa[j]])){
							in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k]),na.rm=T)
						}
					}
					
					if(in_flow >= delta[rem_gene]){
						predict[rem_gene,rem_k]  <- rnorm(1,active_mu[rem_gene],active_sd[rem_gene])
					}
					else{
						predict[rem_gene,rem_k]  <- rnorm(1,inactive_mu[rem_gene],inactive_sd[rem_gene])
					}
				}
			}
		}
	}
	else if (muPgk==T){
	
		for (ent in 1:dim(rem_entries)[1]){
			rem_gene=rem_entries[ent,1]
			rem_k=rem_entries[ent,2]
			
			res=vecMatch(rem_entries[ent,1:2], inact_entries)
			
			# if the removed entry is an inactive node due to some knockdown, then predict as inactivet=rem_entries[ent,3]
			if (any(res)==TRUE){
				predict[rem_gene,rem_k] <- rnorm(1,inactive_mu[rem_gene,rem_k],inactive_sd[rem_gene,rem_k])
			}
			else{
				pa <- which(adja[,rem_gene]!=0)
				# if there are no parents: rem_gene is root node and thus is considered to be active since it is not silenced
				if (length(pa) == 0){
					predict[rem_gene,rem_k]  <- rnorm(1,active_mu[rem_gene,rem_k],active_sd[rem_gene,rem_k]) 
				}
				else{
					in_flow <- baseline[rem_gene]
					
					for(j in 1:length(pa)){
						if (!is.na(obs[pa[j],rem_k]) &  (obs[pa[j],rem_k]>= delta[pa[j],rem_k])){
							in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k]),na.rm=T)
						}
					}
					
					if(in_flow >= delta[rem_gene,rem_k]){
						predict[rem_gene,rem_k]  <- rnorm(1,active_mu[rem_gene,rem_k],active_sd[rem_gene,rem_k])
					}
					else{
						predict[rem_gene,rem_k]  <- rnorm(1,inactive_mu[rem_gene,rem_k],inactive_sd[rem_gene,rem_k])
					}
				}
			}
		}
	}
	
	
	return(predict)
}


calcPredictionKfoldCV_dyn_new <-function(b,n,K,adja,baseline,obs,delta,rem_entries,active_mu,active_sd,inactive_mu,inactive_sd,muPgene=FALSE,muPgk=FALSE,muPgt=FALSE,muPgkt=FALSE)
{
	# activation matrix is the same regardless of time point
	act_mat <- calcActivation(adja,b,n,K)
	inact_entries = which(act_mat==0, arr.ind=T) # returns (i,k)
	predict = obs
	
	if (muPgene==F & muPgk==F & muPgt==F & muPgkt==F){
	
		for (ent in 1:dim(rem_entries)[1]){
			rem_gene=rem_entries[ent,1]
			rem_k=rem_entries[ent,2]
			rem_t=rem_entries[ent,3]
			
			res=vecMatch(rem_entries[ent,1:2], inact_entries)
			
			# if the removed entry is an inactive node due to some knockdown, then predict as inactivet=rem_entries[ent,3]
			if (any(res)==TRUE){
				predict[rem_gene,rem_k,rem_t] <- rnorm(1,inactive_mu,inactive_sd)
			}
			else{
				pa <- which(adja[,rem_gene]!=0)
				# if there are no parents: rem_gene is root node and thus is considered to be active since it is not silenced
				if (length(pa) == 0){
					predict[rem_gene,rem_k,rem_t] <- rnorm(1,active_mu,active_sd) 
				}
				else{
					in_flow <- baseline[rem_gene]
					
					for(j in 1:length(pa)){
						if (!is.na(obs[pa[j],rem_k,rem_t-1]) &  (obs[pa[j],rem_k,rem_t-1]>= delta[pa[j]])){
							in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
						}
					}
					
					if(in_flow >= delta[rem_gene]){
						predict[rem_gene,rem_k,rem_t] <- rnorm(1,active_mu,active_sd)
					}
					else{
						predict[rem_gene,rem_k,rem_t] <- rnorm(1,inactive_mu,inactive_sd)
					}
				}
			}
		}
	}
	
	else if (muPgene==T){
	
		for (ent in 1:dim(rem_entries)[1]){
			rem_gene=rem_entries[ent,1]
			rem_k=rem_entries[ent,2]
			rem_t=rem_entries[ent,3]
			
			res=vecMatch(rem_entries[ent,1:2], inact_entries)
			
			# if the removed entry is an inactive node due to some knockdown, then predict as inactivet=rem_entries[ent,3]
			if (any(res)==TRUE){
				predict[rem_gene,rem_k,rem_t] <- rnorm(1,inactive_mu[rem_gene],inactive_sd[rem_gene])
			}
			else{
				pa <- which(adja[,rem_gene]!=0)
				# if there are no parents: rem_gene is root node and thus is considered to be active since it is not silenced
				if (length(pa) == 0){
					predict[rem_gene,rem_k,rem_t] <- rnorm(1,active_mu[rem_gene],active_sd[rem_gene]) 
				}
				else{
					in_flow <- baseline[rem_gene]
					
					for(j in 1:length(pa)){
						if (!is.na(obs[pa[j],rem_k,rem_t-1]) &  (obs[pa[j],rem_k,rem_t-1]>= delta[pa[j]])){
							in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
						}
					}
					
					if(in_flow >= delta[rem_gene]){
						predict[rem_gene,rem_k,rem_t] <- rnorm(1,active_mu[rem_gene],active_sd[rem_gene])
					}
					else{
						predict[rem_gene,rem_k,rem_t] <- rnorm(1,inactive_mu[rem_gene],inactive_sd[rem_gene])
					}
				}
			}
		}
	}
	else if (muPgk==T){
	
		for (ent in 1:dim(rem_entries)[1]){
			rem_gene=rem_entries[ent,1]
			rem_k=rem_entries[ent,2]
			rem_t=rem_entries[ent,3]
			
			res=vecMatch(rem_entries[ent,1:2], inact_entries)
			
			# if the removed entry is an inactive node due to some knockdown, then predict as inactivet=rem_entries[ent,3]
			if (any(res)==TRUE){
				predict[rem_gene,rem_k,rem_t] <- rnorm(1,inactive_mu[rem_gene,rem_k],inactive_sd[rem_gene,rem_k])
			}
			else{
				pa <- which(adja[,rem_gene]!=0)
				# if there are no parents: rem_gene is root node and thus is considered to be active since it is not silenced
				if (length(pa) == 0){
					predict[rem_gene,rem_k,rem_t] <- rnorm(1,active_mu[rem_gene,rem_k],active_sd[rem_gene,rem_k]) 
				}
				else{
					in_flow <- baseline[rem_gene]
					
					for(j in 1:length(pa)){
						if (!is.na(obs[pa[j],rem_k,rem_t-1]) &  (obs[pa[j],rem_k,rem_t-1]>= delta[pa[j],rem_k])){
							in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
						}
					}
					
					if(in_flow >= delta[rem_gene,rem_k]){
						predict[rem_gene,rem_k,rem_t] <- rnorm(1,active_mu[rem_gene,rem_k],active_sd[rem_gene,rem_k])
					}
					else{
						predict[rem_gene,rem_k,rem_t] <- rnorm(1,inactive_mu[rem_gene,rem_k],inactive_sd[rem_gene,rem_k])
					}
				}
			}
		}
	}
	
	else if (muPgt==T){
	
		for (ent in 1:dim(rem_entries)[1]){
			rem_gene=rem_entries[ent,1]
			rem_k=rem_entries[ent,2]
			rem_t=rem_entries[ent,3]
			
			res=vecMatch(rem_entries[ent,1:2], inact_entries)
			
			# if the removed entry is an inactive node due to some knockdown, then predict as inactivet=rem_entries[ent,3]
			if (any(res)==TRUE){
				predict[rem_gene,rem_k,rem_t] <- rnorm(1,inactive_mu[rem_gene, rem_t],inactive_sd[rem_gene, rem_t])
			}
			else{
				pa <- which(adja[,rem_gene]!=0)
				# if there are no parents: rem_gene is root node and thus is considered to be active since it is not silenced
				if (length(pa) == 0){
					predict[rem_gene,rem_k,rem_t] <- rnorm(1,active_mu[rem_gene, rem_t],active_sd[rem_gene, rem_t]) 
				}
				else{
					in_flow <- baseline[rem_gene]
					
					for(j in 1:length(pa)){
						if (!is.na(obs[pa[j],rem_k,rem_t-1]) &  (obs[pa[j],rem_k,rem_t-1]>= delta[pa[j],rem_t-1])){
							in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
						}
					}
					
					if(in_flow >= delta[rem_gene]){
						predict[rem_gene,rem_k,rem_t] <- rnorm(1,active_mu[rem_gene, rem_t],active_sd[rem_gene, rem_t])
					}
					else{
						predict[rem_gene,rem_k,rem_t] <- rnorm(1,inactive_mu[rem_gene, rem_t],inactive_sd[rem_gene, rem_t])
					}
				}
			}
		}
	}
	
	else if (muPgkt==T){
	
		for (ent in 1:dim(rem_entries)[1]){
			rem_gene=rem_entries[ent,1]
			rem_k=rem_entries[ent,2]
			rem_t=rem_entries[ent,3]
			
			res=vecMatch(rem_entries[ent,1:2], inact_entries)
			
			# if the removed entry is an inactive node due to some knockdown, then predict as inactivet=rem_entries[ent,3]
			if (any(res)==TRUE){
				predict[rem_gene,rem_k,rem_t] <- rnorm(1,inactive_mu[rem_gene, rem_k, rem_t],inactive_sd[rem_gene, rem_k, rem_t])
			}
			else{
				pa <- which(adja[,rem_gene]!=0)
				# if there are no parents: rem_gene is root node and thus is considered to be active since it is not silenced
				if (length(pa) == 0){
					predict[rem_gene,rem_k,rem_t] <- rnorm(1,active_mu[rem_gene, rem_k, rem_t],active_sd[rem_gene, rem_k, rem_t]) 
				}
				else{
					in_flow <- baseline[rem_gene]
					
					for(j in 1:length(pa)){
						if (!is.na(obs[pa[j],rem_k,rem_t-1]) &  (obs[pa[j],rem_k,rem_t-1]>= delta[pa[j],rem_k, rem_t-1])){
							in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
						}
					}
					
					if(in_flow >= delta[rem_gene,rem_k,rem_t]){
						predict[rem_gene,rem_k,rem_t] <- rnorm(1,active_mu[rem_gene, rem_k, rem_t],active_sd[rem_gene, rem_k, rem_t])
					}
					else{
						predict[rem_gene,rem_k,rem_t] <- rnorm(1,inactive_mu[rem_gene, rem_k, rem_t],inactive_sd[rem_gene, rem_k, rem_t])
					}
				}
			}
		}
	}
	
	return(predict)
}

	
vecMatch <- function(vec, mat) {
	out <- apply(mat, 1, function(mat, vec) isTRUE(all.equal(mat, vec)), vec)
	return(out)
}
