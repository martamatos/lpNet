getObsMat <-
function(act_mat,active_mu,active_sd,inactive_mu,inactive_sd,geneState=NULL){
  obs_mat <- act_mat
  
  if (length(geneState!=0)){
		obs_mat[(geneState==1 & act_mat==1)] <- rnorm(sum(geneState==1 & act_mat==1),active_mu,active_sd)
		obs_mat[(geneState==0 | act_mat==0)] <- rnorm(sum(geneState==0 | act_mat==0),inactive_mu,inactive_sd)
	}
	else{
		obs_mat[act_mat==1] <- rnorm(sum(act_mat==1),active_mu,active_sd)
		obs_mat[act_mat==0] <- rnorm(sum(act_mat==0),inactive_mu,inactive_sd)
	}

  return(obs_mat)
}


#getObsMat_ddepn <-
#function(act_mat,active_mu,active_sd,inactive_mu,inactive_sd){
  
#  obs_mat <- act_mat
  
#  act_entries = which(act_mat==1, arr.ind=T)
#  inact_entries = which(act_mat==0, arr.ind=T)
  
#  obs_mat[act_entries] = rnorm(sum(act_mat==1),active_mu[act_entries],active_sd[act_entries])
#  obs_mat[inact_entries] = rnorm(sum(act_mat==0),inactive_mu[inact_entries],inactive_sd[inact_entries])

#  return(obs_mat)
#}

#getObsMat_sahin03 <-
#function(act_mat,active_mu,active_sd,inactive_mu,inactive_sd){
  
#  obs_mat <- act_mat
  
#  act_entries = which(act_mat==1, arr.ind=T)
#  inact_entries = which(act_mat==0, arr.ind=T)
  
#  obs_mat[act_entries] = rnorm(sum(act_mat==1),active_mu[act_entries],active_sd[act_entries])
#  obs_mat[inact_entries] = rnorm(sum(act_mat==0),inactive_mu[inact_entries],inactive_sd[inact_entries])

#  return(obs_mat)
#}


#getObsMat_2lev <-
#function(act_mat,active_mu,active_sd,inactive_mu,inactive_sd, inactivePknock_mu, inactivePknock_sd, geneState=NULL){
#  obs_mat <- act_mat

#  if (length(geneState!=0)){
#		obs_mat[(geneState==1 & act_mat==1)] <- rnorm(sum(geneState==1 & act_mat==1),active_mu,active_sd)
#		obs_mat[(xor(geneState==0, act_mat==0))] <- rnorm(sum(xor(geneState==0, act_mat==0)),inactive_mu,inactive_sd)
#		obs_mat[(geneState==0 & act_mat==0)] <- rnorm(sum( geneState==0 & act_mat==0 ),inactivePknock_mu,inactivePknock_sd)
#	}
#	else{
#		obs_mat[act_mat==1] <- rnorm(sum(act_mat==1),active_mu,active_sd)
#		obs_mat[act_mat==0] <- rnorm(sum(act_mat==0),inactive_mu,inactive_sd)
#	}

#  return(obs_mat)
#}
