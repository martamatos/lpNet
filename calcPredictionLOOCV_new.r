#
# LOOCV prediction for LP original
#

calcPredictionLOOCV_LP <-function(b,n,K,adja,baseline,obs,delta,rem_gene, rem_kd,active_mu,active_sd,inactive_mu,inactive_sd,muPgene=FALSE,muPgk=FALSE)
{
	kds <- matrix(b,nrow=n,ncol=K)
	
	sil_gene_ids <- which(kds[,rem_kd]==0)
	
	if (muPgene==F & muPgk==F){
		# if removed observation is from a gene which has been silenced: prediction comes from inactivation
		if(rem_gene %in% sil_gene_ids){
			predict <- rnorm(1,inactive_mu,inactive_sd)
		}
		# else: in_flow: sum of all parents of rem_gene after knockdown rem_kd times the weights
		else{
			pa <- which(adja[,rem_gene]!=0)
			# if there are no parents: rem_gene is root node and thus is considered to be active since it is not silenced
			if(length(pa)==0){
				predict <- rnorm(1,active_mu,active_sd)
			}
			else{
				# calculate in_flow
				in_flow <- 0
				for(j in 1:length(pa)){
					in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_kd]),na.rm=T)
				}
				# if in_flow >= delta[rem_gene] its active, otherwise not
				if(in_flow>= delta[rem_gene]){
					predict <- rnorm(1,active_mu,active_sd)
				}
				else predict <- rnorm(1,inactive_mu,inactive_sd)
			}
		}
	}
	else if (muPgene==T){
		if(rem_gene %in% sil_gene_ids){
			predict <- rnorm(1,inactive_mu[rem_gene],inactive_sd[rem_gene])
		}
		# else: in_flow: sum of all parents of rem_gene after knockdown rem_kd times the weights
		else{
			pa <- which(adja[,rem_gene]!=0)
			# if there are no parents: rem_gene is root node and thus is considered to be active since it is not silenced
			if(length(pa)==0){
				predict <- rnorm(1,active_mu[rem_gene],active_sd[rem_gene])
			}
			else{
				# calculate in_flow
				in_flow <- 0
				for(j in 1:length(pa)){
					in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_kd]),na.rm=T)
				}
				# if in_flow >= delta[rem_gene] its active, otherwise not
				if(in_flow>= delta[rem_gene]){
					predict <- rnorm(1,active_mu[rem_gene],active_sd[rem_gene])
				}
				else predict <- rnorm(1,inactive_mu[rem_gene],inactive_sd[rem_gene])
			}
		}
	}
	else if (muPgk==T){
		if(rem_gene %in% sil_gene_ids){
			predict <- rnorm(1,inactive_mu[rem_gene,rem_kd],inactive_sd[rem_gene,rem_kd])
		}
		# else: in_flow: sum of all parents of rem_gene after knockdown rem_kd times the weights
		else{
			pa <- which(adja[,rem_gene]!=0)
			# if there are no parents: rem_gene is root node and thus is considered to be active since it is not silenced
			if(length(pa)==0){
				predict <- rnorm(1,active_mu[rem_gene,rem_kd],active_sd[rem_gene,rem_kd])
			}
			else{
				# calculate in_flow
				in_flow <- 0
				for(j in 1:length(pa)){
					in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_kd]),na.rm=T)
				}
				# if in_flow >= delta[rem_gene] its active, otherwise not
				if(in_flow>= delta[rem_gene,rem_kd]){
					predict <- rnorm(1,active_mu[rem_gene,rem_kd],active_sd[rem_gene,rem_kd])
				}
				else predict <- rnorm(1,inactive_mu[rem_gene,rem_kd],inactive_sd[rem_gene,rem_kd])
			}
		}
	
	}
}
	
#
# LOOCV prediction for no discretized model
#
calcPredictionLOOCV_dyn <-function(b,n,K,adja,baseline,obs,delta,rem_gene, rem_k, rem_t,active_mu,active_sd,inactive_mu,inactive_sd,muPgene=FALSE,muPgk=FALSE,muPgt=FALSE,muPgkt=FALSE)
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
				if (is.na(baseline[rem_gene])) base <- 0
				else base <- baseline[rem_gene]
				
				if( base >= delta[rem_gene]){
					predict <- rnorm(1,active_mu,active_sd) 
				}
				else{
					predict <- rnorm(1,inactive_mu,inactive_sd)
				}
			}
			else{
				if (is.na(baseline[rem_gene])) in_flow <- 0
				else in_flow <- baseline[rem_gene]
				
				for(j in 1:length(pa)){
					in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
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
				if (is.na(baseline[rem_gene])) base <- 0
				else base <- baseline[rem_gene]
				
				if (base >= delta[rem_gene]){
					predict <- rnorm(1,active_mu[rem_gene],active_sd[rem_gene]) 
				}
				else{
					predict <- rnorm(1,inactive_mu[rem_gene],inactive_sd[rem_gene])
				}
			}
			else{
				if (is.na(baseline[rem_gene])) in_flow <- 0
				else in_flow <- baseline[rem_gene]
				
				for(j in 1:length(pa)){
					in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
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
				if (is.na(baseline[rem_gene])) base <- 0
				else base <- baseline[rem_gene]
				
				if (base >= delta[rem_gene,rem_k]){
					predict <- rnorm(1,active_mu[rem_gene,rem_k],active_sd[rem_gene,rem_k]) 
				}
				else{
					predict <- rnorm(1,inactive_mu[rem_gene,rem_k],inactive_sd[rem_gene,rem_k])
				}
			}
			else{
				if (is.na(baseline[rem_gene])) in_flow <- 0
				else in_flow <- baseline[rem_gene]
				
				for(j in 1:length(pa)){
					in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
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
				if (is.na(baseline[rem_gene])) base <- 0
				else base <- baseline[rem_gene]
				
				if (base >= delta[rem_gene,rem_t]){
					predict <- rnorm(1,active_mu[rem_gene,rem_t],active_sd[rem_gene,rem_t]) 
				}
				else{
					predict <- rnorm(1,inactive_mu[rem_gene,rem_t],inactive_sd[rem_gene,rem_t])
				}
			}
			else{
				if (is.na(baseline[rem_gene])) in_flow <- 0
				else in_flow <- baseline[rem_gene]
				
				for(j in 1:length(pa)){
					in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
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
				if (is.na(baseline[rem_gene])) base <- 0
				else base <- baseline[rem_gene]
				
				if (base >= delta[rem_gene, rem_k,rem_t]){
					predict <- rnorm(1,active_mu[rem_gene, rem_k,rem_t],active_sd[rem_gene, rem_k,rem_t]) 
				}
				else{
					predict <- rnorm(1,inactive_mu[rem_gene, rem_k,rem_t],inactive_sd[rem_gene, rem_k,rem_t])
				}
			}
			else{
				if (is.na(baseline[rem_gene])) in_flow <- 0
				else in_flow <- baseline[rem_gene]
				
				for(j in 1:length(pa)){
					in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
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

#
# LOOCV prediction for half discretized model
#
	
calcPredictionLOOCV_dyn_disc <-function(b,n,K,adja,baseline,obs,delta,rem_gene, rem_k, rem_t,active_mu,active_sd,inactive_mu,inactive_sd,muPgene=FALSE,muPgk=FALSE,muPgt=FALSE,muPgkt=FALSE)
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
				if (is.na(baseline[rem_gene])) base <- 0
				else base <- baseline[rem_gene]
				
				if( base >= delta[rem_gene]){
					predict <- rnorm(1,active_mu,active_sd) 
				}
				else{
					predict <- rnorm(1,inactive_mu,inactive_sd)
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
				if (is.na(baseline[rem_gene])) base <- 0
				else base <- baseline[rem_gene]
				
				if (base >= delta[rem_gene]){
					predict <- rnorm(1,active_mu[rem_gene],active_sd[rem_gene]) 
				}
				else{
					predict <- rnorm(1,inactive_mu[rem_gene],inactive_sd[rem_gene])
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
				if (is.na(baseline[rem_gene])) base <- 0
				else base <- baseline[rem_gene]
				
				if (base >= delta[rem_gene,rem_k]){
					predict <- rnorm(1,active_mu[rem_gene,rem_k],active_sd[rem_gene,rem_k]) 
				}
				else{
					predict <- rnorm(1,inactive_mu[rem_gene,rem_k],inactive_sd[rem_gene,rem_k])
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
				if (is.na(baseline[rem_gene])) base <- 0
				else base <- baseline[rem_gene]
				
				if (base >= delta[rem_gene,rem_t]){
					predict <- rnorm(1,active_mu[rem_gene,rem_t],active_sd[rem_gene,rem_t]) 
				}
				else{
					predict <- rnorm(1,inactive_mu[rem_gene,rem_t],inactive_sd[rem_gene,rem_t])
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
				if (is.na(baseline[rem_gene])) base <- 0
				else base <- baseline[rem_gene]
				
				if (base >= delta[rem_gene, rem_k,rem_t]){
					predict <- rnorm(1,active_mu[rem_gene, rem_k,rem_t],active_sd[rem_gene, rem_k,rem_t]) 
				}
				else{
					predict <- rnorm(1,inactive_mu[rem_gene, rem_k,rem_t],inactive_sd[rem_gene, rem_k,rem_t])
				}
			}
			else{
				if (is.na(baseline[rem_gene])) in_flow <- 0
				else in_flow <- baseline[rem_gene]
				
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
