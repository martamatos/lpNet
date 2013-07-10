getAdja <-function(result,numnodes)
{
  edges <- cbind(result$solution,names(result$objective)) 
#  print("edges")
#  print(edges)
  exist_edge <- edges[which(edges[,1]!=0),]
  adja <- matrix(0,nrow=numnodes,ncol=numnodes)

  if(!is.null(dim(exist_edge)))
  {
	id_neg <- grep("w-_",exist_edge[,2],fixed=T)
	id_pos <- grep("w+_",exist_edge[,2],fixed=T)
	
		if(length(id_pos)>0)
		{
			pos_edge <- gsub("w+_","",exist_edge[id_pos,2],fixed=T)
			pos_edge <- strsplit(pos_edge,"_")
			for(i in 1:length(id_pos))
			{
				adja[as.integer(pos_edge[[i]][1]),as.integer(pos_edge[[i]][2])] <- as.double(exist_edge[id_pos[i],1])
			}
		}
	
		if(length(id_neg)>0){  
			neg_edge <- gsub("w-_","",exist_edge[id_neg,2],fixed=T)
			neg_edge <- strsplit(neg_edge,"_")
		
			for(i in 1:length(id_neg))
			{
				id1 <- as.integer(neg_edge[[i]][1])
				id2 <- as.integer(neg_edge[[i]][2])
				# check whether there is already a positive edge: if yes: take the one with the higher weight
				if(adja[id1,id2]!=0)
				{
					val_neg <- as.double(exist_edge[id_neg[i],1])
					# if negative is bigger: take negative
					if(adja[id1,id2]<val_neg) adja[id1,id2] <- -val_neg
					# if both are equal: take average: 0
					if(adja[id1,id2]==val_neg) adja[id1,id2] <- 0
					# if positive is bigger: take positive -> thus, do nothing
				}
				else adja[id1,id2] <- -as.double(exist_edge[id_neg[i],1])
			}
		}
  }
  return(adja)
}

