getProbMatrix <- function(edges_all,numnodes,septype="->")
{

	
	nTotal = dim(edges_all)[1]
	nCol = dim(edges_all)[2]
	pPos = rep(NA, nCol)
	pNeg = rep(NA, nCol)
	
	for (j in 1:nCol){
		pPos[j] = sum(edges_all[,j]==1)/nTotal
		pNeg[j] = sum(edges_all[,j]==-1)/nTotal
	}


  sample <- matrix(0,nrow=numnodes,ncol=numnodes)
  colnames(sample) <- rownames(sample) <- seq(1,numnodes)
  edgelist <- strsplit(colnames(edges_all),septype,fixed=T)
	
	
  for(i in 1:nCol)
  {
		id1 <- which(edgelist[[i]][1]==rownames(sample))
		id2 <- which(edgelist[[i]][2]==colnames(sample))

			if( pPos[i] > pNeg[i] ){
				sample[id1,id2] <- pPos[i]
			}
			else if( pPos[i] < pNeg[i] ){
				sample[id1,id2] <- -pNeg[i]
			}
			else{
				sample[id1,id2] <- 0
			}
  }

  return(sample)
}

getProbMatrix02 <- function(edges_all,numnodes,septype="->")
{
#	print("edges alL 1 ")
#	print(edges_all)

#	edges_all2 = getLOOCVActiveEdges(edges_all, numnodes)
#	print("edges all2 ")
#	print(edges_all2)
	
	nTotal = dim(edges_all)[1]
	nCol = dim(edges_all)[2]
	pPos = rep(NA, nCol)
	pNeg = rep(NA, nCol)
	
	for (j in 1:nCol){
		pPos[j] = sum(edges_all[,j]==1)/nTotal
		pNeg[j] = sum(edges_all[,j]==-1)/nTotal
	}
	
#	print("bla1")
#	print(pPos)
#	print(pNeg)
#	print(length(pPos))
#	print(length(pNeg))

  sample <- matrix(0,nrow=numnodes,ncol=numnodes)
  colnames(sample) <- rownames(sample) <- seq(1,numnodes)
  edgelist <- strsplit(colnames(edges_all),septype,fixed=T)
	
	
  for(i in 1:nCol)
  {
		id1 <- which(edgelist[[i]][1]==rownames(sample))
		id2 <- which(edgelist[[i]][2]==colnames(sample))

			if( pPos[i] >= pNeg[i] ){
				sample[id1,id2] <- pPos[i]
			}
			else{
				sample[id1,id2] <- -pNeg[i]
			}
			
  }
#  print(sample)
##  stop("A")
  return(sample)
}

