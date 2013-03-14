getEdgeAnnot <- function(n,allpos=FALSE,T_=NULL)
{
  if(!allpos)
  {
		if ( length(T_)>0 ){	# create annotations for nonIterative model
			ret = c()
			for (t in 1:T_){
				ret <- c(ret, paste("w+",rep(seq(1,n),rep(n,n)),rep(seq(1,n),n),t,sep="_"), # pos variables
								paste("w-",rep(seq(1,n),rep(n,n)),rep(seq(1,n),n),t,sep="_"))  # neg variables
			}
			ret = c(ret, paste("w",seq(1,n),"^",0,sep="_"))
		}
		else{ # create annotations for iterative model
			ret <- c(paste("w+",rep(seq(1,n),rep(n,n)),rep(seq(1,n),n),sep="_"), # pos variables
							paste("w-",rep(seq(1,n),rep(n,n)),rep(seq(1,n),n),sep="_"),  # neg variables
							paste("w",seq(1,n),"^",0,sep="_"))
		}
  }
  else # if all.pos
  {
		if ( length(T_)>0 ){	# create annotations for nonIterative model
			ret = c()
			for (t in 1:T_){
				ret <- c(ret, paste("w+",rep(seq(1,n),rep(n,n)),rep(seq(1,n),n),t,sep="_")) # pos variables
			}
			ret = c(ret, paste("w",seq(1,n),"^",0,sep="_"))
		}
		else{	# create annotations for iterative model
			ret <- c(paste("w+",rep(seq(1,n),rep(n,n)),rep(seq(1,n),n),sep="_"), # pos variables
            paste("w",seq(1,n),"^",0,sep="_"))
    }
  }
  return(ret)
}
