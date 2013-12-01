#
# kfoldCV prediction for LP original
#

calcPredictionKfoldCV_LP <-function(b,n,K,adja,baseline,obs,delta,rem_entries,rem_entries_vec,active_mu,active_sd,inactive_mu,inactive_sd,mu_type=NULL)
{
  # which genes are silenced in removed observation
  act_mat <- calcActivation(adja,b,n,K)
  predict <- getObsMat(act_mat,active_mu,active_sd,inactive_mu,inactive_sd)
  return(predict)
}



#
# kfoldCV prediction for not discretized model
#

calcPredictionKfoldCV_dyn <-function(b,n,K,adja,baseline,obs,delta,rem_entries,rem_entries_vec,active_mu,active_sd,inactive_mu,inactive_sd,mu_type)
{
	# activation matrix is the same regardless of time point
	act_mat <- calcActivation_dyn(adja,b,n,K)
	inact_entries = which(act_mat==0) 
	predict = obs

	if (dim(rem_entries)[1]>0){
		if ( mu_type == "single"){
		
			for (ent in 1:dim(rem_entries)[1]){
				rem_gene=rem_entries[ent,1]
				rem_k=rem_entries[ent,2]
				rem_t=rem_entries[ent,3]
				
				rem_ent_test = rem_entries_vec[ent]%%(n*K)
				if (rem_ent_test == 0) rem_ent_test = n*K # when the removed entry is in the last position the remainder of the division is zero and the test fails
				res= (rem_ent_test %in% inact_entries)
				
			
				# if the removed entry is an inactive node due to some knockdown, then predict as inactivet=rem_entries[ent,3]
				if (res==TRUE){
					predict[rem_gene,rem_k,rem_t] <- rnorm(1,inactive_mu,inactive_sd)
				}
				else{
					pa <- which(adja[,rem_gene]!=0)
					# if there are no parents: rem_gene is root node and thus is considered to be active or inactive according to its baseline value
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
					
						flagNA = 0
						flagBreak = 0
						
						for(j in 1:length(pa))
						{
							if (is.na(obs[pa[j],rem_k,rem_t-1])){
								flagNA = 1
								
								if (adja[pa[j],rem_gene] < 0 & act_mat[pa[j], rem_k]==1){
									predict[rem_gene,rem_k,rem_t] <- NA
									flagBreak = 1
									break 
								}
							}
							in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
						}
						
						if (flagBreak==1) next
						
						if (flagNA == 1 & in_flow < delta[rem_gene]){
							predict[rem_gene,rem_k,rem_t] <- NA
							next
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
		
		else if ( mu_type == "perGene"){
	
			for (ent in 1:dim(rem_entries)[1]){
				rem_gene=rem_entries[ent,1]
				rem_k=rem_entries[ent,2]
				rem_t=rem_entries[ent,3]
				
				rem_ent_test = rem_entries_vec[ent]%%(n*K)
				if (rem_ent_test == 0) rem_ent_test = n*K # when the removed entry is in the last position the remainder of the division is zero and the test fails
				res= (rem_ent_test %in% inact_entries)
				
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
						
						flagNA = 0
						flagBreak = 0
						
						for(j in 1:length(pa))
						{
							if (is.na(obs[pa[j],rem_k,rem_t-1])){
								flagNA = 1
								
								if (adja[pa[j],rem_gene] < 0 & act_mat[pa[j], rem_k]==1){
									predict[rem_gene,rem_k,rem_t] <- NA
									flagBreak = 1
									break 
								}
							}
							in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
						}
						
						if (flagBreak==1) next
						
						if (flagNA == 1 & in_flow < delta[rem_gene]){
							predict[rem_gene,rem_k,rem_t] <- NA
							next
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
		else if (mu_type == "perGeneExp"){
		
			for (ent in 1:dim(rem_entries)[1]){
				rem_gene=rem_entries[ent,1]
				rem_k=rem_entries[ent,2]
				rem_t=rem_entries[ent,3]
				
				rem_ent_test = rem_entries_vec[ent]%%(n*K)
				if (rem_ent_test == 0) rem_ent_test = n*K # when the removed entry is in the last position the remainder of the division is zero and the test fails
				res= (rem_ent_test %in% inact_entries)
				
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
						
						flagNA = 0
						flagBreak = 0
						
						for(j in 1:length(pa))
						{
							if (is.na(obs[pa[j],rem_k,rem_t-1])){
								flagNA = 1
								
								if (adja[pa[j],rem_gene] < 0 & act_mat[pa[j], rem_k]==1){
									predict[rem_gene,rem_k,rem_t] <- NA
									flagBreak = 1
									break 
								}
							}
							in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
						}
						
						if (flagBreak==1) next
						
						if (flagNA == 1 & in_flow < delta[rem_gene,rem_k]){
							predict[rem_gene,rem_k,rem_t] <- NA
							next
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
		
		else if (mu_type == "perGeneTime"){
		
			for (ent in 1:dim(rem_entries)[1]){
				rem_gene=rem_entries[ent,1]
				rem_k=rem_entries[ent,2]
				rem_t=rem_entries[ent,3]
				
				rem_ent_test = rem_entries_vec[ent]%%(n*K)
				if (rem_ent_test == 0) rem_ent_test = n*K # when the removed entry is in the last position the remainder of the division is zero and the test fails
				res= (rem_ent_test %in% inact_entries)
				
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
												
						flagNA = 0
						flagBreak = 0
						
						for(j in 1:length(pa))
						{
							if (is.na(obs[pa[j],rem_k,rem_t-1])){
								flagNA = 1
								
								if (adja[pa[j],rem_gene] < 0 & act_mat[pa[j], rem_k]==1){
									predict[rem_gene,rem_k,rem_t] <- NA
									flagBreak = 1
									break 
								}
							}
							in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
						}
						
						if (flagBreak==1) next
						
						if (flagNA == 1 & in_flow < delta[rem_gene,rem_t]){
							predict[rem_gene,rem_k,rem_t] <- NA
							next
						}
						
						if(in_flow >= delta[rem_gene,rem_t]){
							predict[rem_gene,rem_k,rem_t] <- rnorm(1,active_mu[rem_gene,rem_t],active_sd[rem_gene,rem_t]) 
						}
						else{
							predict[rem_gene,rem_k,rem_t] <- rnorm(1,inactive_mu[rem_gene,rem_t],inactive_sd[rem_gene,rem_t])
						}
					}
				}
			}
		}
		
		else if (mu_type == "perGeneExpTime"){
		
			for (ent in 1:dim(rem_entries)[1]){
				rem_gene=rem_entries[ent,1]
				rem_k=rem_entries[ent,2]
				rem_t=rem_entries[ent,3]
				
				rem_ent_test = rem_entries_vec[ent]%%(n*K)
				if (rem_ent_test == 0) rem_ent_test = n*K # when the removed entry is in the last position the remainder of the division is zero and the test fails
				res= (rem_ent_test %in% inact_entries)
				
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
						
						flagNA = 0
						flagBreak = 0
						
						for(j in 1:length(pa))
						{
							if (is.na(obs[pa[j],rem_k,rem_t-1])){
								flagNA = 1
								
								if (adja[pa[j],rem_gene] < 0 & act_mat[pa[j], rem_k]==1){
									predict[rem_gene,rem_k,rem_t] <- NA
									flagBreak = 1
									break 
								}
							}
							in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
						}
						
						if (flagBreak==1) next
						
						if (flagNA == 1 & in_flow < delta[rem_gene,rem_k,rem_t]){
							predict[rem_gene,rem_k,rem_t] <- NA
							next
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
calcPredictionKfoldCV_dyn_disc <-function(b,n,K,adja,baseline,obs,delta,rem_entries,rem_entries_vec,active_mu,active_sd,inactive_mu,inactive_sd,mu_type)
{
	# activation matrix is the same regardless of time point
	act_mat <- calcActivation_dyn(adja,b,n,K)
	inact_entries = which(act_mat==0) 
	predict = obs

	if (dim(rem_entries)[1]>0){
		if ( mu_type == "single"){
		
			for (ent in 1:dim(rem_entries)[1]){
				rem_gene=rem_entries[ent,1]
				rem_k=rem_entries[ent,2]
				rem_t=rem_entries[ent,3]
				
				rem_ent_test = rem_entries_vec[ent]%%(n*K)
				if (rem_ent_test == 0) rem_ent_test = n*K # when the removed entry is in the last position the remainder of the division is zero and the test fails
				res= (rem_ent_test %in% inact_entries)
				
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
						
						flagNA = 0
						flagBreak = 0
						
						for(j in 1:length(pa))
						{
							if (is.na(obs[pa[j],rem_k,rem_t-1])){
								flagNA = 1
								
								if (adja[pa[j],rem_gene] < 0 & act_mat[pa[j], rem_k]==1){
									predict[rem_gene,rem_k,rem_t] <- NA
									flagBreak = 1
									break 
								}
							}
							else if (obs[pa[j],rem_k,rem_t-1]>= delta[pa[j]]){
								in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
							}
						}
				
						if (flagBreak==1) next
						
						if (flagNA == 1 & in_flow < delta[rem_gene]){
							predict[rem_gene,rem_k,rem_t] <- NA
							next
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
		
		else if (mu_type == "perGene"){

			for (ent in 1:dim(rem_entries)[1]){
				rem_gene=rem_entries[ent,1]
				rem_k=rem_entries[ent,2]
				rem_t=rem_entries[ent,3]
				
				rem_ent_test = rem_entries_vec[ent]%%(n*K)
				if (rem_ent_test == 0) rem_ent_test = n*K # when the removed entry is in the last position the remainder of the division is zero and the test fails
				res= (rem_ent_test %in% inact_entries)
				
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
						
						flagNA = 0
						flagBreak = 0
						
						for(j in 1:length(pa))
						{
							if (is.na(obs[pa[j],rem_k,rem_t-1])){
								flagNA = 1
								
								if (adja[pa[j],rem_gene] < 0 & act_mat[pa[j], rem_k]==1){
									predict[rem_gene,rem_k,rem_t] <- NA
									flagBreak = 1
									break 
								}
							}
							else if (obs[pa[j],rem_k,rem_t-1]>= delta[pa[j]]){
								in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
							}
						}
						
						if (flagBreak==1) next
						
						if (flagNA == 1 & in_flow < delta[rem_gene]){
							predict[rem_gene,rem_k,rem_t] <- NA
							next
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
		else if (mu_type == "perGeneExp"){
		
			for (ent in 1:dim(rem_entries)[1]){
				rem_gene=rem_entries[ent,1]
				rem_k=rem_entries[ent,2]
				rem_t=rem_entries[ent,3]
				
				rem_ent_test = rem_entries_vec[ent]%%(n*K)
				if (rem_ent_test == 0) rem_ent_test = n*K # when the removed entry is in the last position the remainder of the division is zero and the test fails
				res= (rem_ent_test %in% inact_entries)
				
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
						
						flagNA = 0
						flagBreak = 0
						
						for(j in 1:length(pa))
						{
							if (is.na(obs[pa[j],rem_k,rem_t-1])){
								flagNA = 1
								
								if (adja[pa[j],rem_gene] < 0 & act_mat[pa[j], rem_k]==1){
									predict[rem_gene,rem_k,rem_t] <- NA
									flagBreak = 1
									break 
								}
							}
							else if (obs[pa[j],rem_k,rem_t-1]>= delta[pa[j],rem_k]){
								in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
							}
						}
						
						if (flagBreak==1) next
						
						if (flagNA == 1 & in_flow < delta[rem_gene,rem_k]){
							predict[rem_gene,rem_k,rem_t] <- NA
							next
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
		
		else if (mu_type == "perGeneTime"){
		
			for (ent in 1:dim(rem_entries)[1]){
				rem_gene=rem_entries[ent,1]
				rem_k=rem_entries[ent,2]
				rem_t=rem_entries[ent,3]
				
				rem_ent_test = rem_entries_vec[ent]%%(n*K)
				if (rem_ent_test == 0) rem_ent_test = n*K # when the removed entry is in the last position the remainder of the division is zero and the test fails
				res= (rem_ent_test %in% inact_entries)
				
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
						
						flagNA = 0
						flagBreak = 0
						
						for(j in 1:length(pa))
						{
							if (is.na(obs[pa[j],rem_k,rem_t-1])){
								flagNA = 1
								
								if (adja[pa[j],rem_gene] < 0 & act_mat[pa[j], rem_k]==1){
									predict[rem_gene,rem_k,rem_t] <- NA
									flagBreak = 1
									break 
								}
							}
							else if (obs[pa[j],rem_k,rem_t-1]>= delta[pa[j],rem_t-1]){
								in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
							}
						}
						
						if (flagBreak==1) next
						
						if (flagNA == 1 & in_flow < delta[rem_gene,rem_t]){
							predict[rem_gene,rem_k,rem_t] <- NA
							next
						}			
						
						if(in_flow >= delta[rem_gene,rem_t]){
							predict[rem_gene,rem_k,rem_t] <- rnorm(1,active_mu[rem_gene,rem_t],active_sd[rem_gene,rem_t])
						}
						else{
							predict[rem_gene,rem_k,rem_t] <- rnorm(1,inactive_mu[rem_gene,rem_t],inactive_sd[rem_gene,rem_t])
						}
					}
				}
			}
		}
		
		else if (mu_type == "perGeneExpTime"){
		
			for (ent in 1:dim(rem_entries)[1]){
				rem_gene=rem_entries[ent,1]
				rem_k=rem_entries[ent,2]
				rem_t=rem_entries[ent,3]
				
				rem_ent_test = rem_entries_vec[ent]%%(n*K)
				if (rem_ent_test == 0) rem_ent_test = n*K # when the removed entry is in the last position the remainder of the division is zero and the test fails
				res= (rem_ent_test %in% inact_entries)
				
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
						
						flagNA = 0
						flagBreak = 0
						
						for(j in 1:length(pa))
						{
							if (is.na(obs[pa[j],rem_k,rem_t-1])){
								flagNA = 1
								
								if (adja[pa[j],rem_gene] < 0 & act_mat[pa[j], rem_k]==1){
									predict[rem_gene,rem_k,rem_t] <- NA
									flagBreak = 1
									break 
								}
							}
							else if (obs[pa[j],rem_k,rem_t-1]>= delta[pa[j],rem_k,rem_t-1]){
								in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
							}
						}
						
						if (flagBreak==1) next
						
						if (flagNA == 1 & in_flow < delta[rem_gene,rem_k,rem_t]){
							predict[rem_gene,rem_k,rem_t] <- NA
							next
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
calcPredictionKfoldCV_dyn_disc_dream8 <-function(b,n,K,adja,baseline,obs,delta,rem_entries,rem_entries_vec,active_mu,active_sd,inactive_mu,inactive_sd,mu_type)
{
	# activation matrix is the same regardless of time point
	act_mat <- calcActivation_dyn(adja,b,n,K)
	inact_entries = which(act_mat==0) 
	predict = obs

	if (dim(rem_entries)[1]>0){
		if ( mu_type == "single"){
		
			for (ent in 1:dim(rem_entries)[1]){
				rem_gene=rem_entries[ent,1]
				rem_k=rem_entries[ent,2]
				rem_t=rem_entries[ent,3]
				
				rem_ent_test = rem_entries_vec[ent]%%(n*K)
				if (rem_ent_test == 0) rem_ent_test = n*K # when the removed entry is in the last position the remainder of the division is zero and the test fails
				res= (rem_ent_test %in% inact_entries)
				
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
						
						flagNA = 0
						flagBreak = 0
						
						for(j in 1:length(pa))
						{
							if (is.na(obs[pa[j],rem_k,rem_t-1])){
								flagNA = 1
								
								if (adja[pa[j],rem_gene] < 0 & act_mat[pa[j], rem_k]==1){
									predict[rem_gene,rem_k,rem_t] <- NA
									flagBreak = 1
									break 
								}
							}
							else if (obs[pa[j],rem_k,rem_t-1]>= delta[pa[j]] & b[(rem_k-1)*n + pa[j]]==1){
								in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
							}
						}
				
						if (flagBreak==1) next
						
						if (flagNA == 1 & in_flow < delta[rem_gene]){
							predict[rem_gene,rem_k,rem_t] <- NA
							next
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
		
		else if (mu_type == "perGene"){

			for (ent in 1:dim(rem_entries)[1]){
				rem_gene=rem_entries[ent,1]
				rem_k=rem_entries[ent,2]
				rem_t=rem_entries[ent,3]
				
				rem_ent_test = rem_entries_vec[ent]%%(n*K)
				if (rem_ent_test == 0) rem_ent_test = n*K # when the removed entry is in the last position the remainder of the division is zero and the test fails
				res= (rem_ent_test %in% inact_entries)
				
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
						
						flagNA = 0
						flagBreak = 0
						
						for(j in 1:length(pa))
						{
							if (is.na(obs[pa[j],rem_k,rem_t-1])){
								flagNA = 1
								
								if (adja[pa[j],rem_gene] < 0 & act_mat[pa[j], rem_k]==1){
									predict[rem_gene,rem_k,rem_t] <- NA
									flagBreak = 1
									break 
								}
							}
							else if (obs[pa[j],rem_k,rem_t-1]>= delta[pa[j]] & b[(rem_k-1)*n + pa[j]]==1){
								in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
							}
						}
						
						if (flagBreak==1) next
						
						if (flagNA == 1 & in_flow < delta[rem_gene]){
							predict[rem_gene,rem_k,rem_t] <- NA
							next
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
		else if (mu_type == "perGeneExp"){
		
			for (ent in 1:dim(rem_entries)[1]){
				rem_gene=rem_entries[ent,1]
				rem_k=rem_entries[ent,2]
				rem_t=rem_entries[ent,3]
				
				rem_ent_test = rem_entries_vec[ent]%%(n*K)
				if (rem_ent_test == 0) rem_ent_test = n*K # when the removed entry is in the last position the remainder of the division is zero and the test fails
				res= (rem_ent_test %in% inact_entries)
				
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
						
						flagNA = 0
						flagBreak = 0
						
						for(j in 1:length(pa))
						{
							if (is.na(obs[pa[j],rem_k,rem_t-1])){
								flagNA = 1
								
								if (adja[pa[j],rem_gene] < 0 & act_mat[pa[j], rem_k]==1){
									predict[rem_gene,rem_k,rem_t] <- NA
									flagBreak = 1
									break 
								}
							}
							else if (obs[pa[j],rem_k,rem_t-1]>= delta[pa[j],rem_k] & b[(rem_k-1)*n + pa[j]]==1){
								in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
							}
						}
						
						if (flagBreak==1) next
						
						if (flagNA == 1 & in_flow < delta[rem_gene,rem_k]){
							predict[rem_gene,rem_k,rem_t] <- NA
							next
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
		
		else if (mu_type == "perGeneTime"){
		
			for (ent in 1:dim(rem_entries)[1]){
				rem_gene=rem_entries[ent,1]
				rem_k=rem_entries[ent,2]
				rem_t=rem_entries[ent,3]
				
				rem_ent_test = rem_entries_vec[ent]%%(n*K)
				if (rem_ent_test == 0) rem_ent_test = n*K # when the removed entry is in the last position the remainder of the division is zero and the test fails
				res= (rem_ent_test %in% inact_entries)
				
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
						
						flagNA = 0
						flagBreak = 0
						
						for(j in 1:length(pa))
						{
							if (is.na(obs[pa[j],rem_k,rem_t-1])){
								flagNA = 1
								
								if (adja[pa[j],rem_gene] < 0 & act_mat[pa[j], rem_k]==1){
									predict[rem_gene,rem_k,rem_t] <- NA
									flagBreak = 1
									break 
								}
							}
							else if (obs[pa[j],rem_k,rem_t-1]>= delta[pa[j],rem_t-1] & b[(rem_k-1)*n + pa[j]]==1){
								in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
							}
						}
						
						if (flagBreak==1) next
						
						if (flagNA == 1 & in_flow < delta[rem_gene,rem_t]){
							predict[rem_gene,rem_k,rem_t] <- NA
							next
						}			
						
						if(in_flow >= delta[rem_gene,rem_t]){
							predict[rem_gene,rem_k,rem_t] <- rnorm(1,active_mu[rem_gene,rem_t],active_sd[rem_gene,rem_t])
						}
						else{
							predict[rem_gene,rem_k,rem_t] <- rnorm(1,inactive_mu[rem_gene,rem_t],inactive_sd[rem_gene,rem_t])
						}
					}
				}
			}
		}
		
		else if (mu_type == "perGeneExpTime"){
		
			for (ent in 1:dim(rem_entries)[1]){
				rem_gene=rem_entries[ent,1]
				rem_k=rem_entries[ent,2]
				rem_t=rem_entries[ent,3]
				
				rem_ent_test = rem_entries_vec[ent]%%(n*K)
				if (rem_ent_test == 0) rem_ent_test = n*K # when the removed entry is in the last position the remainder of the division is zero and the test fails
				res= (rem_ent_test %in% inact_entries)
				
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
						
						flagNA = 0
						flagBreak = 0
						
						for(j in 1:length(pa))
						{
							if (is.na(obs[pa[j],rem_k,rem_t-1])){
								flagNA = 1
								
								if (adja[pa[j],rem_gene] < 0 & act_mat[pa[j], rem_k]==1){
									predict[rem_gene,rem_k,rem_t] <- NA
									flagBreak = 1
									break 
								}
							}
							else if (obs[pa[j],rem_k,rem_t-1]>= delta[pa[j],rem_k,rem_t-1] & b[(rem_k-1)*n + pa[j]]==1){
								in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_k,rem_t-1]),na.rm=T)
							}
						}
						
						if (flagBreak==1) next
						
						if (flagNA == 1 & in_flow < delta[rem_gene,rem_k,rem_t]){
							predict[rem_gene,rem_k,rem_t] <- NA
							next
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
