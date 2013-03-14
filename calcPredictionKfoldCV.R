calcPredictionKfoldCV <-
function(adja,b,n,K,active_mu,active_sd,inactive_mu,inactive_sd){
  # which genes are silenced in removed observation
  act_mat <- calcActivation(adja,b,n,K)
  predict <- getObsMat(act_mat,active_mu,active_sd,inactive_mu,inactive_sd)
  return(predict)
}

calcPredictionKfoldCV_dyn <-
function(T_,adja,b,n,K,active_mu,active_sd,inactive_mu,inactive_sd){
  # which genes are silenced in removed observation
  
#  print("calc act")
#  print(adja)
#  print(b)
#  print(n)
#  print(K)
  act_mat <- calcActivation(adja,b,n,K)
#  print("act_mat")
#  print(act_mat)
  predict = array(NA, c(n,K,T_))
  for (t in 1:T_){
#		print("adja[,,t]")
#		print(adja[,,t])
#		print(geneState[,,t])
		
		predict[,,t] <- getObsMat(act_mat,active_mu,active_sd,inactive_mu,inactive_sd)
	}
#	print("predict")
#	print(predict)
  return(predict)
}


calcPredictionKfoldCV_nonIterative <-
function(T_,adja,b,n,K,active_mu,active_sd,inactive_mu,inactive_sd){
  # which genes are silenced in removed observation
  
  predict = array(NA, c(n,K,T_))
  for (t in 1:T_){
#		print("adja[,,t]")
#		print(adja[,,t])
#		print(geneState[,,t])
		act_mat <- calcActivation(adja[,,t],b,n,K)
		predict[,,t] <- getObsMat(act_mat,active_mu,active_sd,inactive_mu,inactive_sd)
	}
#	print("predict")
#	print(predict)
  return(predict)
}
