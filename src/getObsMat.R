#
# generate observation matrix for simulated data
#
getObsMat <-
function(act_mat,active_mu,active_sd,inactive_mu,inactive_sd,geneState=NULL){
  obs_mat <- act_mat
    
  if (!is.null(geneState)){
		obs_mat[geneState==1] <- rnorm(sum(geneState==1),active_mu,active_sd)
		obs_mat[geneState==0] <- rnorm(sum(geneState==0),inactive_mu,inactive_sd)
	}
	else{
		obs_mat[act_mat==1] <- rnorm(sum(act_mat==1),active_mu,active_sd)
		obs_mat[act_mat==0] <- rnorm(sum(act_mat==0),inactive_mu,inactive_sd)
	}

  return(obs_mat)
}


# to test different data ranges for each gene
getObsMat_range <-
function(act_mat,active_mu,active_sd,inactive_mu,inactive_sd,geneState=NULL){
  
  obs_mat <- geneState
  
  print(geneState)
  print(dim(geneState))
  for (gene in 1:length(geneState)){
		if (geneState[gene] == 1){
			obs_mat[gene] <- rnorm(1,active_mu[gene],active_sd)
		}
		else if (geneState[gene] == 0){
			obs_mat[gene] <- rnorm(1,inactive_mu[gene],inactive_sd)
		}
	}

  return(obs_mat)
}
