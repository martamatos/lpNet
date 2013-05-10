calcPredictionLOOCV_dyn <-function(rem_t,kds,adja,obs,delta,rem_kd,rem_gene,active_mu,active_sd,inactive_mu,inactive_sd)
{
  # which genes are silenced in removed observation
  sil_gene_ids <- which(kds[,rem_kd,rem_t]==0)

  # if removed observation is from a gene which has been silenced: prediction comes from inactivation
  if(rem_gene %in% sil_gene_ids){
#		predict <- rnorm(1,inactive_mu,inactive_sd)
		if (length(inactive_mu) > 1){
#			print("len 1")
			predict <- rnorm(1,inactive_mu[rem_gene],inactive_sd[rem_gene])
#			print(rem_gene)
#			print(inactive_mu)
#			print(inactive_sd)
		}
		else{
			predict <- rnorm(1,inactive_mu,inactive_sd)
		}
	}
	
	# else: in_flow: sum of all parents of rem_gene after knockdown rem_kd times the weights
	else{
		pa <- which(adja[,rem_gene]!=0)
		# if there are no parents: rem_gene is root node and thus is considered to be active since it is not silenced
		if(length(pa)==0){
#			predict <- rnorm(1,active_mu,active_sd)
			if (length(active_mu) > 1){
#				print("len 1")
				predict <- rnorm(1,active_mu[rem_gene],active_sd[rem_gene])
#				print(rem_gene)
#				print(active_mu)
#				print(active_sd)
			}
			else{
#				print("len >1")
				predict <- rnorm(1,active_mu,active_sd)
			}
		}
		else{
			# calculate in_flow
			in_flow <- 0
			for(j in 1:length(pa)){
				in_flow <- sum(in_flow,(adja[pa[j],rem_gene]*obs[pa[j],rem_kd,rem_t]),na.rm=T)
			}
			# if in_flow >= delta[rem_gene] its active, otherwise not
			if(in_flow >= delta[rem_gene]){
				if (length(active_mu) > 1){
#					print("len 1")
					predict <- rnorm(1,active_mu[rem_gene],active_sd[rem_gene])
#					print(rem_gene)
#					print(active_mu)
#					print(active_sd)
				}
				else{
#					print("len >1")
					predict <- rnorm(1,active_mu,active_sd)
				}
			}
			else{
				if (length(inactive_mu) > 1){
#					print("len 1")
					predict <- rnorm(1,inactive_mu[rem_gene],inactive_sd[rem_gene])
#					print(rem_gene)
#					print(inactive_mu)
#					print(inactive_sd)
				}
				else{
#					print("len >1")
					predict <- rnorm(1,inactive_mu,inactive_sd)
				}
			}
			
#			 predict <- rnorm(1,inactive_mu,inactive_sd)
		}
	}
  return(predict)
}
