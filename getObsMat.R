getObsMat <-
function(act_mat,active_mu,active_sd,inactive_mu,inactive_sd,geneState=NULL){
  obs_mat <- act_mat
#  print("act_mat")
#  print(act_mat)
#  print(geneState)
##  print(obs_mat)
##  print(geneState)
#  print(active_mu)
#  print(active_sd)
#  print(inactive_mu)
#  print(inactive_sd)
  
  if (length(geneState!=0)){
		obs_mat[(geneState==1 & act_mat==1)] <- rnorm(sum(geneState==1 & act_mat==1),active_mu,active_sd)
		obs_mat[(geneState==0 | act_mat==0)] <- rnorm(sum(geneState==0 | act_mat==0),inactive_mu,inactive_sd)
	}
	else{
		obs_mat[act_mat==1] <- rnorm(sum(act_mat==1),active_mu,active_sd)
		obs_mat[act_mat==0] <- rnorm(sum(act_mat==0),inactive_mu,inactive_sd)
	}
#  print("obs_mat")
#  print(obs_mat)

  return(obs_mat)
}
