calcRangeLambda <-
function(delta,obs,stepsize){
  l <- 0
  for(k in 1:length(delta)){
		l <- sum(c(l,obs[k,]<delta[k]),na.rm=T)
  }
  if(l==0) l <- 1
  tmp <- round(l*var(c(obs),na.rm=T),digits=2)
  if(tmp<stepsize) tmp <- stepsize
  return(lambda=c(seq(0,tmp,by=stepsize)))
}

calcRangeLambda_nonIterative <-
function(delta,obs,stepsize){

#	print("cacl lambda")
#	print(obs)
#	print(delta)
  l <- 0
  for(k in 1:length(delta)){
		l <- sum(c(l,obs[k,,]<delta[k]),na.rm=T)
  }
  if(l==0) l <- 1
  tmp <- round(l*var(c(obs),na.rm=T),digits=2)
  if(tmp<stepsize) tmp <- stepsize
  return(lambda=c(seq(0,tmp,by=stepsize)))
}

calcRangeLambda_nonIterative_ddepn <-
function(delta,obs,stepsize){

#	print("cacl lambda")
#	print(obs)
#	print(delta)
  l <- 0
  for(i in 1:dim(delta)[1]){
		for(k in 1:dim(delta)[2]){
			l <- sum(c(l,obs[i,k,]<delta[i,k]),na.rm=T)
		}
  }
  if(l==0) l <- 1
  tmp <- round(l*var(c(obs),na.rm=T),digits=2)
  if(tmp<stepsize) tmp <- stepsize
  return(lambda=c(seq(0,tmp,by=stepsize)))
}

calcRangeLambda_nonIterative_sahin03 <-
function(delta,obs,stepsize){

#	print("cacl lambda")
#	print(obs)
#	print(delta)
  l <- 0
  for(i in 1:dim(delta)[1]){
		for(k in 1:dim(delta)[2]){
			l <- sum(c(l,obs[i,,k]<delta[i,k]),na.rm=T)
		}
  }
  if(l==0) l <- 1
  tmp <- round(l*var(c(obs),na.rm=T),digits=2)
  if(tmp<stepsize) tmp <- stepsize
  return(lambda=c(seq(0,tmp,by=stepsize)))
}
