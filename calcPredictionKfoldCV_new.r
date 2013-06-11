#
# kfoldCV prediction for LP original
#

calcPredictionKfoldCV_LP <-function(b,n,K,adja,baseline,obs,delta,rem_entries,rem_entries_vec,active_mu,active_sd,inactive_mu,inactive_sd,muPgene=FALSE,muPgk=FALSE)
{
  # which genes are silenced in removed observation
  act_mat <- calcActivation(adja,b,n,K)
  predict <- getObsMat(act_mat,active_mu,active_sd,inactive_mu,inactive_sd)
  return(predict)
}



#
# kfoldCV prediction for not discretized model
#

calcPredictionKfoldCV_dyn <-function(b,n,K,adja,baseline,obs,delta,rem_entries,rem_entries_vec,active_mu,active_sd,inactive_mu,inactive_sd,muPgene=FALSE,muPgk=FALSE,muPgt=FALSE,muPgkt=FALSE)
{
	# activation matrix is the same regardless of time point
	act_mat <- calcActivation(adja,b,n,K)
	inact_entries = which(act_mat==0) 
	predict = obs

	if (dim(rem_entries)[1]>0){
		if (muPgene==F & muPgk==F & muPgt==F & muPgkt==F){
		
			for (ent in 1:dim(rem_entries)[1]){
				rem_gene=rem_entries[ent,1]
				rem_k=rem_entries[ent,2]
				rem_t=rem_entries[ent,3]
				
				res= (rem_entries_vec[ent] %in% inact_entries)
				
				# if the removed entry is an inactive node due to some knockdown, then predict as inactivet=rem_entries[ent,3]
				if (res==TRUE){
					predict[rem_gene,rem_k,rem_t] <- rnorm(1,inactive_mu,inactive_sd)
				}
				else{
					pa <- which(adja[,rem_gene]!=0)
					# if there are no parents: rem_gene is root node and thus is considered to be active since it is not silenced
					if (length(pa) == 0){
						if (is.na(baseline[rem_gene])) base <- 0
						else base <- baseline[rem_gene]
						
						if (base >= delta[rem_gene]){
							predict[rem_gene,rem_k,rem_t] <- rnorm(1,active_mu,active_sd) 
						}
						else{
							predict[rem_gene,rem_k,rem_t] <- rnorm(1,inactive_mu,inactive_sd)
						}
					}
					else{
						if (is.na(baseline[rem_gene])) in_flow <- 0
						else in_flow <- baseline[rem_gene]
						
						for(j in 1:length(pa)){
							in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
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
				
				res= (rem_entries_vec[ent] %in% inact_entries)
				
				# if the removed entry is an inactive node due to some knockdown, then predict as inactivet=rem_entries[ent,3]
				if (res==TRUE){
					predict[rem_gene,rem_k,rem_t] <- rnorm(1,inactive_mu[rem_gene],inactive_sd[rem_gene])
				}
				else{
					pa <- which(adja[,rem_gene]!=0)
					# if there are no parents: rem_gene is root node and thus is considered to be active since it is not silenced
					if (length(pa) == 0){
						if (is.na(baseline[rem_gene])) base <- 0
						else base <- baseline[rem_gene]
						
						if (base >= delta[rem_gene]){
							predict[rem_gene,rem_k,rem_t] <- rnorm(1,active_mu[rem_gene],active_sd[rem_gene]) 
						}
						else{
							predict[rem_gene,rem_k,rem_t] <- rnorm(1,inactive_mu[rem_gene],inactive_sd[rem_gene])
						}
					}
					else{
						
						if (is.na(baseline[rem_gene])) in_flow <- 0
						else in_flow <- baseline[rem_gene]
						
						for(j in 1:length(pa)){
							in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
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
				
				res= (rem_entries_vec[ent] %in% inact_entries)
				
				# if the removed entry is an inactive node due to some knockdown, then predict as inactivet=rem_entries[ent,3]
				if (res==TRUE){
					predict[rem_gene,rem_k,rem_t] <- rnorm(1,inactive_mu[rem_gene,rem_k],inactive_sd[rem_gene,rem_k])
				}
				else{
					pa <- which(adja[,rem_gene]!=0)
					# if there are no parents: rem_gene is root node and thus is considered to be active since it is not silenced
					
					if (length(pa) == 0){
						if (is.na(baseline[rem_gene])) base <- 0
						else base <- baseline[rem_gene]
						
						if (base >= delta[rem_gene,rem_k]){
							predict[rem_gene,rem_k,rem_t] <- rnorm(1,active_mu[rem_gene,rem_k],active_sd[rem_gene,rem_k]) 
						}
						else{
							predict[rem_gene,rem_k,rem_t] <- rnorm(1,inactive_mu[rem_gene,rem_k],inactive_sd[rem_gene,rem_k])
						}
					}
					else{
						if (is.na(baseline[rem_gene])) in_flow <- 0
						else in_flow <- baseline[rem_gene]
						for(j in 1:length(pa)){
							in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
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
				
				res= (rem_entries_vec[ent] %in% inact_entries)
				
				# if the removed entry is an inactive node due to some knockdown, then predict as inactivet=rem_entries[ent,3]
				if (res==TRUE){
					predict[rem_gene,rem_k,rem_t] <- rnorm(1,inactive_mu[rem_gene, rem_t],inactive_sd[rem_gene, rem_t])
				}
				else{
					pa <- which(adja[,rem_gene]!=0)
					# if there are no parents: rem_gene is root node and thus is considered to be active since it is not silenced
					if (length(pa) == 0){
						if (is.na(baseline[rem_gene])) base <- 0
						else base <- baseline[rem_gene]
						
						if (base >= delta[rem_gene, rem_t]){
							predict[rem_gene,rem_k,rem_t] <- rnorm(1,active_mu[rem_gene, rem_t],active_sd[rem_gene, rem_t]) 
						}
						else{
							predict[rem_gene,rem_k,rem_t] <- rnorm(1,inactive_mu[rem_gene, rem_t],inactive_sd[rem_gene, rem_t])
						}
					}
					else{
						if (is.na(baseline[rem_gene])) in_flow <- 0
						else in_flow <- baseline[rem_gene]
						
						for(j in 1:length(pa)){
							in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
						}
						
						if(in_flow >= delta[rem_gene, rem_t]){
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
				
				res= (rem_entries_vec[ent] %in% inact_entries)
				
				# if the removed entry is an inactive node due to some knockdown, then predict as inactivet=rem_entries[ent,3]
				if (res==TRUE){
					predict[rem_gene,rem_k,rem_t] <- rnorm(1,inactive_mu[rem_gene, rem_k, rem_t],inactive_sd[rem_gene, rem_k, rem_t])
				}
				else{
					pa <- which(adja[,rem_gene]!=0)
					# if there are no parents: rem_gene is root node and thus is considered to be active since it is not silenced
					if (length(pa) == 0){
						if (is.na(baseline[rem_gene])) base <- 0
						else base <- baseline[rem_gene]
						
						if (base >= delta[rem_gene,rem_k,rem_t]){
							predict[rem_gene,rem_k,rem_t] <- rnorm(1,active_mu[rem_gene, rem_k, rem_t],active_sd[rem_gene, rem_k, rem_t]) 
						}
						else{
							predict[rem_gene,rem_k,rem_t] <- rnorm(1,inactive_mu[rem_gene, rem_k, rem_t],inactive_sd[rem_gene, rem_k, rem_t])
						}
					}
					else{
						if (is.na(baseline[rem_gene])) in_flow <- 0
						else in_flow <- baseline[rem_gene]
						
						for(j in 1:length(pa)){
							in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
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
		
	}
	return(predict)
}



#
# kfoldCV prediction for half discretized model
#
calcPredictionKfoldCV_dyn_disc <-function(b,n,K,adja,baseline,obs,delta,rem_entries,rem_entries_vec,active_mu,active_sd,inactive_mu,inactive_sd,muPgene=FALSE,muPgk=FALSE,muPgt=FALSE,muPgkt=FALSE)
{
	# activation matrix is the same regardless of time point
	act_mat <- calcActivation(adja,b,n,K)
	inact_entries = which(act_mat==0) 
	predict = obs

	if (dim(rem_entries)[1]>0){
		if (muPgene==F & muPgk==F & muPgt==F & muPgkt==F){
		
			for (ent in 1:dim(rem_entries)[1]){
				rem_gene=rem_entries[ent,1]
				rem_k=rem_entries[ent,2]
				rem_t=rem_entries[ent,3]
				
				res= (rem_entries_vec[ent] %in% inact_entries)
				
				# if the removed entry is an inactive node due to some knockdown, then predict as inactivet=rem_entries[ent,3]
				if (res==TRUE){
					predict[rem_gene,rem_k,rem_t] <- rnorm(1,inactive_mu,inactive_sd)
				}
				else{
					pa <- which(adja[,rem_gene]!=0)
					# if there are no parents: rem_gene is root node and thus is considered to be active since it is not silenced
					if (length(pa) == 0){
						if (is.na(baseline[rem_gene])) base <- 0
						else base <- baseline[rem_gene]
						
						if (base >= delta[rem_gene]){
							predict[rem_gene,rem_k,rem_t] <- rnorm(1,active_mu,active_sd) 
						}
						else{
							predict[rem_gene,rem_k,rem_t] <- rnorm(1,inactive_mu,inactive_sd)
						}
					}
					else{
						if (is.na(baseline[rem_gene])) in_flow <- 0
						else in_flow <- baseline[rem_gene]
						
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
				
				res= (rem_entries_vec[ent] %in% inact_entries)
				
				# if the removed entry is an inactive node due to some knockdown, then predict as inactivet=rem_entries[ent,3]
				if (res==TRUE){
					predict[rem_gene,rem_k,rem_t] <- rnorm(1,inactive_mu[rem_gene],inactive_sd[rem_gene])
				}
				else{
					pa <- which(adja[,rem_gene]!=0)
					# if there are no parents: rem_gene is root node and thus is considered to be active since it is not silenced
					if (length(pa) == 0){
						if (is.na(baseline[rem_gene])) base <- 0
						else base <- baseline[rem_gene]
						
						if (base >= delta[rem_gene]){
							predict[rem_gene,rem_k,rem_t] <- rnorm(1,active_mu[rem_gene],active_sd[rem_gene]) 
						}
						else{
							predict[rem_gene,rem_k,rem_t] <- rnorm(1,inactive_mu[rem_gene],inactive_sd[rem_gene])
						}
					}
					else{
						
						if (is.na(baseline[rem_gene])) in_flow <- 0
						else in_flow <- baseline[rem_gene]
						
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
				
				res= (rem_entries_vec[ent] %in% inact_entries)
				
				# if the removed entry is an inactive node due to some knockdown, then predict as inactivet=rem_entries[ent,3]
				if (res==TRUE){
					predict[rem_gene,rem_k,rem_t] <- rnorm(1,inactive_mu[rem_gene,rem_k],inactive_sd[rem_gene,rem_k])
				}
				else{
					pa <- which(adja[,rem_gene]!=0)
					# if there are no parents: rem_gene is root node and thus is considered to be active since it is not silenced
					
					if (length(pa) == 0){
						if (is.na(baseline[rem_gene])) base <- 0
						else base <- baseline[rem_gene]
						
						if (base >= delta[rem_gene,rem_k]){
							predict[rem_gene,rem_k,rem_t] <- rnorm(1,active_mu[rem_gene,rem_k],active_sd[rem_gene,rem_k]) 
						}
						else{
							predict[rem_gene,rem_k,rem_t] <- rnorm(1,inactive_mu[rem_gene,rem_k],inactive_sd[rem_gene,rem_k])
						}
					}
					else{
						if (is.na(baseline[rem_gene])) in_flow <- 0
						else in_flow <- baseline[rem_gene]
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
				
				res= (rem_entries_vec[ent] %in% inact_entries)
				
				# if the removed entry is an inactive node due to some knockdown, then predict as inactivet=rem_entries[ent,3]
				if (res==TRUE){
					predict[rem_gene,rem_k,rem_t] <- rnorm(1,inactive_mu[rem_gene, rem_t],inactive_sd[rem_gene, rem_t])
				}
				else{
					pa <- which(adja[,rem_gene]!=0)
					# if there are no parents: rem_gene is root node and thus is considered to be active since it is not silenced
					if (length(pa) == 0){
						if (is.na(baseline[rem_gene])) base <- 0
						else base <- baseline[rem_gene]
						
						if (base >= delta[rem_gene, rem_t]){
							predict[rem_gene,rem_k,rem_t] <- rnorm(1,active_mu[rem_gene, rem_t],active_sd[rem_gene, rem_t]) 
						}
						else{
							predict[rem_gene,rem_k,rem_t] <- rnorm(1,inactive_mu[rem_gene, rem_t],inactive_sd[rem_gene, rem_t])
						}
					}
					else{
						if (is.na(baseline[rem_gene])) in_flow <- 0
						else in_flow <- baseline[rem_gene]
						
						for(j in 1:length(pa)){
							if (!is.na(obs[pa[j],rem_k,rem_t-1]) &  (obs[pa[j],rem_k,rem_t-1]>= delta[pa[j],rem_t-1])){
								in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
							}
						}
						
						if(in_flow >= delta[rem_gene, rem_t]){
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
				
				res= (rem_entries_vec[ent] %in% inact_entries)
				
				# if the removed entry is an inactive node due to some knockdown, then predict as inactivet=rem_entries[ent,3]
				if (res==TRUE){
					predict[rem_gene,rem_k,rem_t] <- rnorm(1,inactive_mu[rem_gene, rem_k, rem_t],inactive_sd[rem_gene, rem_k, rem_t])
				}
				else{
					pa <- which(adja[,rem_gene]!=0)
					# if there are no parents: rem_gene is root node and thus is considered to be active since it is not silenced
					if (length(pa) == 0){
						if (is.na(baseline[rem_gene])) base <- 0
						else base <- baseline[rem_gene]
						
						if (base >= delta[rem_gene,rem_k,rem_t]){
							predict[rem_gene,rem_k,rem_t] <- rnorm(1,active_mu[rem_gene, rem_k, rem_t],active_sd[rem_gene, rem_k, rem_t]) 
						}
						else{
							predict[rem_gene,rem_k,rem_t] <- rnorm(1,inactive_mu[rem_gene, rem_k, rem_t],inactive_sd[rem_gene, rem_k, rem_t])
						}
					}
					else{
						if (is.na(baseline[rem_gene])) in_flow <- 0
						else in_flow <- baseline[rem_gene]
						
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
		
	}
	return(predict)
}