#
# get the edges annotations
#
getEdgeAnnot <- function(n,allpos=FALSE){
  if(!allpos){
		ret <- c(paste("w+",rep(seq(1,n),rep(n,n)),rep(seq(1,n),n),sep="_"), # pos variables
							paste("w-",rep(seq(1,n),rep(n,n)),rep(seq(1,n),n),sep="_"),  # neg variables
							paste("w",seq(1,n),"^",0,sep="_"))
   }
  else{
   ret <- c(paste("w+",rep(seq(1,n),rep(n,n)),rep(seq(1,n),n),sep="_"), # pos variables
            paste("w",seq(1,n),"^",0,sep="_"))
  }
  return(ret)
}
