calcRangeLambda_LP <- function(delta,obs,stepsize,deltaPk=FALSE){
  l <- 0
 
  if (deltaPk == F){
		for(i in 1:length(delta)){
			l <- sum(c(l,obs[i,]<delta[i]),na.rm=T)
		}
	}
	else if (deltaPk == T){
		for(i in 1:dim(delta)[1]){
			for(k in 1:dim(delta)[2]){
				l <- sum(c(l,obs[i,k]<delta[i,k]),na.rm=T)
			}
		}
	}
	
  if(l==0) l <- 1
  tmp <- round(l*var(c(obs),na.rm=T),digits=2)
  if(tmp<stepsize) tmp <- stepsize
  return(lambda=c(seq(0,tmp,by=stepsize)))
}


calcRangeLambda_dyn <- function(delta,obs,stepsize,deltaPk=FALSE,deltaPt=FALSE,deltaPkt=FALSE){

  l <- 0
  
  if (deltaPk == F & deltaPt==F & deltaPkt==F){
		for(i in 1:length(delta)){
			l <- sum(c(l,obs[i,,]<delta[i]),na.rm=T)
		}
	}
	else if (deltaPk == T){
		for(i in 1:dim(delta)[1]){
			for(k in 1:dim(delta)[2]){
				l <- sum(c(l,obs[i,k,]<delta[i,k]),na.rm=T)
			}
		}
	}
	else if (deltaPt == T){
		for(i in 1:dim(delta)[1]){
			for(t in 1:dim(delta)[2]){
				l <- sum(c(l,obs[i,,t]<delta[i,t]),na.rm=T)
				
			}
		}
	}
	else if (deltaPkt == T){
		for(i in 1:dim(delta)[1]){
			for (k in 1:dim(delta)[2]){
				for(t in 1:dim(delta)[3]){
					l <- sum(c(l,obs[i,k,t]<delta[i,k,t]),na.rm=T)
				}
			}
		}
	}
	
  if(l==0) l <- 1
  tmp <- round(l*var(c(obs),na.rm=T),digits=2)
  if(tmp<stepsize) tmp <- stepsize
  return(lambda=c(seq(0,tmp,by=stepsize)))
}

