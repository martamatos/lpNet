calcPredictionLOOCV_new <-function(b,n,K,adja,baseline,obs,delta,rem_gene, rem_k,active_mu,active_sd,inactive_mu,inactive_sd,muPgene=FALSE,muPgk=FALSE)
{
	# activation matrix is the same regardless of time point
	act_mat <- calcActivation(adja,b,n,K)
	inact_entries = which(act_mat==0, arr.ind=T) # returns (i,k)
	
	res=vecMatch(c(rem_gene, rem_k), inact_entries)
	
	
	
	if (muPgene==F & muPgk==F & muPgt==F & muPgkt==F){
		# if the removed entry is an inactive node due to some knockdown, then predict as inactive
		if (any(res)==TRUE){
			predict <- rnorm(1,inactive_mu,inactive_sd)
		}
		else{
			pa <- which(adja[,rem_gene]!=0)
			# if there are no parents: rem_gene is root node and thus is considered to be active since it is not silenced
			if (length(pa) == 0){
				predict <- rnorm(1,active_mu,active_sd) 
			}
			else{
				in_flow <- baseline[rem_gene]
				
				for(j in 1:length(pa)){
					if (!is.na(obs[pa[j],rem_k]) & (obs[pa[j],rem_k] >= delta[pa[j]])){
						in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k]),na.rm=T)
					}
				}
				
				if(in_flow >= delta[rem_gene]){
					predict <- rnorm(1,active_mu,active_sd)
				}
				else{
					predict <- rnorm(1,inactive_mu,inactive_sd)
				}
			}
		}
	}
	# if there is an in/active_mu/sd per gene
	else if(muPgene==T){
		if (any(res)==TRUE){
			predict <- rnorm(1,inactive_mu[rem_gene],inactive_sd[rem_gene])
		}
		else{
			pa <- which(adja[,rem_gene]!=0)
			# if there are no parents: rem_gene is root node and thus is considered to be active since it is not silenced
			if (length(pa) == 0){
				predict <- rnorm(1,active_mu[rem_gene],active_sd[rem_gene]) 
			}
			else{
				in_flow <- baseline[rem_gene]
				for(j in 1:length(pa)){
					if (!is.na(obs[pa[j],rem_k]) & (obs[pa[j],rem_k] >= delta[pa[j]])){
						in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k]),na.rm=T)
					}
				}
				if(in_flow >= delta[rem_gene]){
					predict <- rnorm(1,active_mu[rem_gene],active_sd[rem_gene])
				}
				else{
					predict <- rnorm(1,inactive_mu[rem_gene],inactive_sd[rem_gene])
				}
			}
		}
	}
	# if there is an in/active_mu/sd and delta per gene per knockdown exp
	else if(muPgk==T){
		# if the removed entry is an inactive node due to some knockdown, then predict as inactive
		if (any(res)==TRUE){
			predict <- rnorm(1,inactive_mu[rem_gene,rem_k],inactive_sd[rem_gene,rem_k])
		}
		else{
			pa <- which(adja[,rem_gene]!=0)
			# if there are no parents: rem_gene is root node and thus is considered to be active since it is not silenced
			if (length(pa) == 0){
				predict <- rnorm(1,active_mu[rem_gene,rem_k],active_sd[rem_gene,rem_k]) 
			}
			else{
				in_flow <- baseline[rem_gene]
				for(j in 1:length(pa)){
					if (!is.na(obs[pa[j],rem_k]) & (obs[pa[j],rem_k] >= delta[pa[j],rem_k])){
						in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k]),na.rm=T)
					}
				}
				if(in_flow >= delta[rem_gene,rem_k]){
					predict <- rnorm(1,active_mu[rem_gene,rem_k],active_sd[rem_gene,rem_k])
				}
				else{
					predict <- rnorm(1,inactive_mu[rem_gene,rem_k],inactive_sd[rem_gene,rem_k])
				}
			}
		}
	}
	
	return(predict)
}
	
	
calcPredictionLOOCV_dyn_new <-function(b,n,K,adja,baseline,obs,delta,rem_gene, rem_k, rem_t,active_mu,active_sd,inactive_mu,inactive_sd,muPgene=FALSE,muPgk=FALSE,muPgt=FALSE,muPgkt=FALSE)
{
	# activation matrix is the same regardless of time point
	act_mat <- calcActivation(adja,b,n,K)
	inact_entries = which(act_mat==0, arr.ind=T) # returns (i,k)
	
	res=vecMatch(c(rem_gene, rem_k), inact_entries)
	
	
	if (muPgene==F & muPgk==F & muPgt==F & muPgkt==F){
		# if the removed entry is an inactive node due to some knockdown, then predict as inactive
		if (any(res)==TRUE){
			predict <- rnorm(1,inactive_mu,inactive_sd)
		}
		else{
			pa <- which(adja[,rem_gene]!=0)
			# if there are no parents: rem_gene is root node and thus is considered to be active since it is not silenced
			if (length(pa) == 0){
				predict <- rnorm(1,active_mu,active_sd) 
			}
			else{
				in_flow <- baseline[rem_gene]
				
				for(j in 1:length(pa)){
					if (!is.na(obs[pa[j],rem_k,rem_t-1]) &  (obs[pa[j],rem_k,rem_t-1]>= delta[pa[j]])){
						in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
					}
				}
				
				if(in_flow >= delta[rem_gene]){
					predict <- rnorm(1,active_mu,active_sd)
				}
				else{
					predict <- rnorm(1,inactive_mu,inactive_sd)
				}
			}
		}
	}
	# if there is an in/active_mu/sd per gene
	else if(muPgene==T){
		if (any(res)==TRUE){
			predict <- rnorm(1,inactive_mu[rem_gene],inactive_sd[rem_gene])
		}
		else{
			pa <- which(adja[,rem_gene]!=0)
			# if there are no parents: rem_gene is root node and thus is considered to be active since it is not silenced
			if (length(pa) == 0){
				predict <- rnorm(1,active_mu[rem_gene],active_sd[rem_gene]) 
			}
			else{
				in_flow <- baseline[rem_gene]
				for(j in 1:length(pa)){
					if (!is.na(obs[pa[j],rem_k,rem_t-1]) &  (obs[pa[j],rem_k,rem_t-1]>= delta[pa[j]])){
						in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
					}
				}
				if(in_flow >= delta[rem_gene]){
					predict <- rnorm(1,active_mu[rem_gene],active_sd[rem_gene])
				}
				else{
					predict <- rnorm(1,inactive_mu[rem_gene],inactive_sd[rem_gene])
				}
			}
		}
	}
	# if there is an in/active_mu/sd and delta per gene per knockdown exp
	else if(muPgk==T){
		# if the removed entry is an inactive node due to some knockdown, then predict as inactive
		if (any(res)==TRUE){
			predict <- rnorm(1,inactive_mu[rem_gene,rem_k],inactive_sd[rem_gene,rem_k])
		}
		else{
			pa <- which(adja[,rem_gene]!=0)
			# if there are no parents: rem_gene is root node and thus is considered to be active since it is not silenced
			if (length(pa) == 0){
				predict <- rnorm(1,active_mu[rem_gene,rem_k],active_sd[rem_gene,rem_k]) 
			}
			else{
				in_flow <- baseline[rem_gene]
				for(j in 1:length(pa)){
					if (!is.na(obs[pa[j],rem_k,rem_t-1]) &  (obs[pa[j],rem_k,rem_t-1]>= delta[pa[j],rem_k])){
						in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
					}
				}
				if(in_flow >= delta[rem_gene,rem_k]){
					predict <- rnorm(1,active_mu[rem_gene,rem_k],active_sd[rem_gene,rem_k])
				}
				else{
					predict <- rnorm(1,inactive_mu[rem_gene,rem_k],inactive_sd[rem_gene,rem_k])
				}
			}
		}
	}
	# if there is an in/active_mu/sd and delta per gene per time point
	else if(muPgt==T){
		# if the removed entry is an inactive node due to some knockdown, then predict as inactive
		if (any(res)==TRUE){
			predict <- rnorm(1,inactive_mu[rem_gene,rem_t],inactive_sd[rem_gene,rem_t])
		}
		else{
			pa <- which(adja[,rem_gene]!=0)
			# if there are no parents: rem_gene is root node and thus is considered to be active since it is not silenced
			if (length(pa) == 0){
				predict <- rnorm(1,active_mu[rem_gene,rem_t],active_sd[rem_gene,rem_t]) 
			}
			else{
				in_flow <- baseline[rem_gene]
				for(j in 1:length(pa)){
					if (!is.na(obs[pa[j],rem_k,rem_t-1]) &  (obs[pa[j],rem_k,rem_t-1]>= delta[pa[j],rem_t-1])){
						in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
					}
				}
				if(in_flow >= delta[rem_gene,rem_t]){
					predict <- rnorm(1,active_mu[rem_gene,rem_t],active_sd[rem_gene,rem_t])
				}
				else{
					predict <- rnorm(1,inactive_mu[rem_gene,rem_t],inactive_sd[rem_gene,rem_t])
				}
			}
		}
	}
	# if there is an in/active_mu/sd and delta per gene per knockdown exp per time point
	else if (muPgkt == T){
		# if the removed entry is an inactive node due to some knockdown, then predict as inactive
		if (any(res)==TRUE){
			predict <- rnorm(1,inactive_mu[rem_gene, rem_k,rem_t],inactive_sd[rem_gene, rem_k,rem_t])
		}
		else{
			pa <- which(adja[,rem_gene]!=0)
			# if there are no parents: rem_gene is root node and thus is considered to be active since it is not silenced
			if (length(pa) == 0){
				predict <- rnorm(1,active_mu[rem_gene, rem_k,rem_t],active_sd[rem_gene, rem_k,rem_t]) 
			}
			else{
				in_flow <- baseline[rem_gene]
				for(j in 1:length(pa)){
					if (!is.na(obs[pa[j],rem_k,rem_t-1]) &  (obs[pa[j],rem_k,rem_t-1]>= delta[pa[j], rem_k,rem_t-1])){
						in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
					}
				}
				if(in_flow >= delta[rem_gene, rem_k,rem_t]){
					predict <- rnorm(1,active_mu[rem_gene, rem_k,rem_t],active_sd[rem_gene, rem_k,rem_t])
				}
				else{
					predict <- rnorm(1,inactive_mu[rem_gene, rem_k,rem_t],inactive_sd[rem_gene, rem_k,rem_t])
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
