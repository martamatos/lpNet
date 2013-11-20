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



calcRangeLambda_dyn_dream8 <- function(delta,obs,stepsize,deltaPk=FALSE,deltaPt=FALSE,deltaPkt=FALSE){

  l <- 0
  
  obst = apply(obs, c(1,2,3), mean)
  
  if (deltaPk == F & deltaPt==F & deltaPkt==F){
		for(i in 1:length(delta)){
			l <- sum(c(l,obst[i,,]<delta[i]),na.rm=T)
		}
	}
	else if (deltaPk == T){
		for(i in 1:dim(delta)[1]){
			for(k in 1:dim(delta)[2]){
				l <- sum(c(l,obst[i,k,]<delta[i,k]),na.rm=T)
			}
		}
	}
	else if (deltaPt == T){
		for(i in 1:dim(delta)[1]){
			for(t in 1:dim(delta)[2]){
				l <- sum(c(l,obst[i,,t]<delta[i,t]),na.rm=T)
				
			}
		}
	}
	else if (deltaPkt == T){
		for(i in 1:dim(delta)[1]){
			for (k in 1:dim(delta)[2]){
				for(t in 1:dim(delta)[3]){
					l <- sum(c(l,obst[i,k,t]<delta[i,k,t]),na.rm=T)
				}
			}
		}
	}
	
  if(l==0) l <- 1
  tmp <- round(l*var(c(obs),na.rm=T),digits=2)
  
  if(tmp<stepsize){
		tmp <- stepsize
		lambda = c(seq(0,tmp,by=stepsize))
	}
	else{
		
		lambda = c()
		
		if (tmp > 0.01){
			lambda = c(seq(0,0.1,by=0.01))
		}
		else{
			lambda = c(seq(0,tmp,by=0.01))
			return(unique(c(lambda,tmp)))
		}
		
		if (tmp > 1){
			lambda = c(lambda, seq(0.1,1,by=0.02))
		}
		else{
			lambda = c(lambda, seq(0.1,tmp,by=0.02))
			return(unique(c(lambda,tmp)))
		}
		
		if (tmp > 2){
			lambda = c(lambda, seq(1,2,by=0.05))
		}
		else{
			lambda = c(lambda, seq(1,tmp,by=0.05))
			return(unique(c(lambda,tmp)))
		}
		
		if (tmp > 10){
			lambda = c(lambda, seq(2,10,by=1))
		}
		else{
			lambda = c(lambda, seq(2,tmp,by=1))
			rreturn(unique(c(lambda,tmp)))
		}
		
		if (tmp > 100){
			lambda = c(lambda, seq(10,100,by=10))
		}
		else{
			lambda = c(lambda, seq(10,tmp,by=10))
			return(unique(c(lambda,tmp)))
		}
		
		if (tmp > 1000){
			lambda = c(lambda, seq(100,1000,by=100))
		}
		else{
			lambda = c(lambda, seq(100,tmp,by=100))
			return(unique(c(lambda,tmp)))
		}
		
		if (tmp > 10000){
			lambda = c(lambda, seq(1000,10000,by=1000))
		}
		else{
			lambda = c(lambda, seq(1000,tmp,by=1000))
			return(unique(c(lambda,tmp)))
		}
	
		while (lambda[length(lambda)] < tmp){
			val = lambda[length(lambda)]*2
			lambda = c(lambda, val)
		}
		
		return(unique(c(lambda,tmp)))
	
	}
  return(lambda)
}


