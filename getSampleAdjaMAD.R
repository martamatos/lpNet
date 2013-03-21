getSampleAdjaMAD <- function(edges_all,numnodes,annot_node,method=median,method2=mad,septype="->")
{
#	print("edges alL 1 ")
#	print(edges_all)

#	edges_all2 = getLOOCVActiveEdges(edges_all, numnodes)
	edges_all2=edges_all

#	print("edges all2 ")
#	print(edges_all2)
	
  edge_med <- apply(edges_all2,2,method,na.rm=T)
  edge_mad <- apply(edges_all2,2,method2,na.rm=T)
  
  sample <- matrix(0,nrow=numnodes,ncol=numnodes)
  colnames(sample) <- rownames(sample) <- seq(1,numnodes)
  edgelist <- strsplit(colnames(edges_all),septype,fixed=T)
	
#	print("edge_med")
#	print(length(edge_med))
#	print(edge_med)
	
	
  for(i in 1:length(edge_med))
  {
		id1 <- which(edgelist[[i]][1]==rownames(sample))
		id2 <- which(edgelist[[i]][2]==colnames(sample))
			
#			print(paste("i",i, sep=" "))
#			print(paste("edge_med[i]",edge_med[i], sep=" "))
#			print(paste("edge_mad[i]",edge_mad[i], sep=" "))
			
			if(abs(edge_med[i])>abs(edge_mad[i])){
				sample[id1,id2] <- edge_med[i]
			}
			
  }
  return(sample)
}




#
# goes through all the edges returned from loocv and substitutes the edges
# in which one of the nodes was removed for NA. thus, these edges won't count
# for the median/mean calculation
#
getLOOCVActiveEdges = function(edges_all, numnodes)
{
	edges_all2 = matrix(0, nrow = numnodes, ncol=numnodes*(numnodes-1))
	rownames(edges_all2) = rownames(edges_all)
	colnames(edges_all2) = colnames(edges_all)
	
	for (i in 1:numnodes)
	{
		for (block in 1:numnodes)
		{
			for (j in 1:numnodes)
			{
				if (j > block) jj = (block-1)*(numnodes-1) + (j-1)
				if (j < block) jj = (block-1)*(numnodes-1) + j
				
				if( ((block == i) || (j == i)) && (j != block)){
#					print(paste("block: ", block))
#					print(paste("i: ", i))
#					print(paste("j: ", j))
#					print(paste("jj: ", jj))
					edges_all2[i,jj] = NA
				}
				else{
					if (j!=block) edges_all2[i,jj] = edges_all[i,jj]
				}
			}
		}
	}
	
	
	return(edges_all2)
}


getSampleAdjaMAD_nonIterative <- function(edges_all,numnodes, annot_node,method=median,method2=mad,septype="->")
{
#	print("edges alL 1 ")
#	print(edges_all)
	
	edges_all2 = edges_all
#	for (t in 1:T_){
#		edges_all2 = rbind(edges_all2, getLOOCVActiveEdges(edges_all, numnodes))
#	}
	

#	print("edges all2 ")
#	print(edges_all2)
	
  edge_med <- apply(edges_all2,2,method,na.rm=T)
  edge_mad <- apply(edges_all2,2,method2,na.rm=T)
  
  sample <- matrix(0,nrow=numnodes,ncol=numnodes)
  colnames(sample) <- rownames(sample) <- annot_node
  edgelist <- strsplit(colnames(edges_all),septype,fixed=T)
	
#	print("edge_med")
#	print(length(edge_med))
#	print(edge_med)
	
	
  for(i in 1:length(edge_med))
  {
		id1 <- which(edgelist[[i]][1]==rownames(sample))
		id2 <- which(edgelist[[i]][2]==colnames(sample))
		
#			print(paste("i",i, sep=" "))
#			print(paste("edge_med[i]",edge_med[i], sep=" "))
#			print(paste("edge_mad[i]",edge_mad[i], sep=" "))
		
		if(abs(edge_med[i])>abs(edge_mad[i])){
			sample[id1,id2] <- edge_med[i]
		}
		
  }
  
  return(sample)
}


