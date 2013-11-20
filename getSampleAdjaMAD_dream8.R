#
# sum up all inferred networks into a single one using method.
#  all edges are included, and the final edge value will be the number
#  of networks in which it was inferred divided by the total number of
#  inferred networks
#
getSampleAdjaMAD_dream8 <- function(edges_all,numnodes,annot_node,method=median,method2=mad,method2Times=1,septype="->")
{

	edges_all2=edges_all
	
	totalEdges = dim(edges_all)[1]
	
  edge_med <- apply(edges_all2,2,method,na.rm=T)
  edge_mad <- apply(edges_all2,2,method2,na.rm=T)
  
  sample <- matrix(0,nrow=numnodes,ncol=numnodes)
  colnames(sample) <- rownames(sample) <- seq(1,numnodes)
  edgelist <- strsplit(colnames(edges_all),septype,fixed=T)

  for(i in 1:length(edgelist))
  {
		id1 <- which(edgelist[[i]][1]==rownames(sample))
		id2 <- which(edgelist[[i]][2]==colnames(sample))


#		if(abs(edge_med[i])>abs(edge_mad[i])){	
		
			posWeight = length(which(edges_all2[,i]>0))/totalEdges
			negWeight = length(which(edges_all2[,i]<0))/totalEdges
		
			if (posWeight > negWeight){
				sample[id1,id2] <- posWeight
			}
			else if (posWeight < negWeight){
				sample[id1,id2] <- -negWeight
			}
#		}
  }
  return(sample)
}


