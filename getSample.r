getSampleMatrix <- function(edges_all,numnodes, method=median,septype="->")
{
#	print("dim getSample matrix")
#  print(edges_all)
  
	annot_node <- seq(1,numnodes)
	edge_med <- apply(edges_all,2,method,na.rm=T)
  sample <- matrix(0,nrow=numnodes,ncol=numnodes)
  colnames(sample) <- rownames(sample) <- annot_node
  edgelist <- strsplit(colnames(edges_all),septype,fixed=T)
  
  for(i in 1:length(edge_med)){
		id1 <- which(edgelist[[i]][1]==rownames(sample))
		id2 <- which(edgelist[[i]][2]==colnames(sample))
		sample[id1,id2] <- edge_med[i]
  }
  
  return(sample)
}
