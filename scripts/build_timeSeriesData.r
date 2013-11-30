getNTimePoints = function(nw_und,nodeNames,n,K)
{
	
	activeNW = list()
	activeNW_temp = matrix(0, nrow=n, ncol=n)
	activeNW[[1]] = activeNW_temp
	activeNW[[2]] = activeNW_temp
	
	geneState = list()
	activeGVec = rep(0,n)
	geneState[[1]] = matrix(rep(activeGVec, K), nrow=n, ncol=K)
	
	in_deg <- apply(abs(nw_und),2,sum)
	parentNodes = which(in_deg==0)
	rootNodes = parentNodes
	nonRootNodes = which(!(seq(1,n) %in% rootNodes))
	
	print(paste("rootNodes ", rootNodes))
	
	
	activeNW_temp = matrix(0, nrow=n, ncol=n)
	
	activeGVec = rep(0,n)
	activeGVec[rootNodes] = 1
	geneState[[2]] = matrix(rep(activeGVec, K), nrow=n, ncol=K)
	
	count = 1
	visitedNodes = vector()
	edgesAdded = vector()
	nEdges = length(which(nw_und!=0))
	
	# the loop ends when all edges have been found or when there are no more children
	while ((length(unique(edgesAdded))<nEdges)){
		
		allChildren = c()

		for (i in parentNodes){
			children = c()
			possChild =which(nw_und[i,]!=0)

			# condition to avoid loops in net construction
			for (j in possChild){
				if (activeNW[[1+count]][i,j]==0){
					children <- c(children, j)
				}
			}
			allChildren = c(allChildren, children)
			activeNW_temp[i, children] = nw_und[i,children]
			edgesAdded = c(edgesAdded, sprintf("%s->%s", i, children))
		}
		
		if (length(allChildren) == 0){
			break
		}

		activeNW[[2+count]] = activeNW_temp
		activeGVec[allChildren] = 1
	
	
		for (i in nonRootNodes){	
			if ( sum(activeNW_temp[,i]) > 0) activeGVec[i] = 1
			else activeGVec[i] = 0						
		}
		for (i in nonRootNodes){
			if (activeGVec[i] == 0) activeNW_temp[i,]=0
		}

		geneState[[2+count]] = matrix(rep(activeGVec, K), nrow=n, ncol=K)

		parentNodes = which(activeGVec > 0)
		count = count + 1
		
		edgesAdded = unique(edgesAdded)

	}
	
#	print(geneState)
#	print(activeNW)
	return(length(geneState))
}

build_timeSeries = function(nw_und,nodeNames,b,n,K)
{

	T_ = getNTimePoints(nw_und,nodeNames,n,1)
	
	act_mat = calcActivation_dyn(nw_und,b,n,K)				
	
	activeNW = list()
	activeNW_temp = matrix(0, nrow=n, ncol=n)
	activeNW[[1]] = activeNW_temp
	activeNW[[2]] = activeNW_temp
	
	geneState = list()
	activeGVec = rep(0,n)
	geneState[[1]] = matrix(rep(activeGVec,K), nrow=n, ncol=K)
	
	for ( t in 2:T_) geneState[[t]] = matrix(NA, nrow=n, ncol=K)
	
	
	for (k in 1:K)
	{
		in_deg <- apply(abs(nw_und),2,sum)
		parentNodes = which(in_deg==0)
		silencedParents = which(parentNodes %in% which(act_mat[,k]==0))
	
		if (length(silencedParents) >= 1) parentNodes = parentNodes[-which(parentNodes %in% which(act_mat[,k]==0))] # to be substituted by -silencedParents
	
		rootNodes = parentNodes
		nonRootNodes = which(!(seq(1,n) %in% rootNodes))
			
		activeNW_temp = matrix(0, nrow=n, ncol=n)
		
		activeGVec = rep(0,n)
		activeGVec[rootNodes] = 1
		activeGVec[which(act_mat[,k]==0)] = 0 # redundant
		geneState[[2]][,k] = activeGVec
		
		count = 1
		visitedNodes = vector()
		edgesAdded = vector()
		nEdges = length(which(nw_und!=0))
		

#		while ((length(unique(edgesAdded))<nEdges)){
		for (bla in 3:T_){
		
			allChildren = c()

			for (i in parentNodes){
				children = c()
				possChild =which(nw_und[i,]!=0)

				# condition to avoid loops in net construction
				for (j in possChild){
					if (activeNW[[1+count]][i,j]==0){
						children <- c(children, j)
					}
				}
				allChildren = c(allChildren, children)
				activeNW_temp[i, children] = nw_und[i,children]
				edgesAdded = c(edgesAdded, sprintf("%s->%s", i, children)) # no need for this anymore
			}
			
#			if (length(allChildren) == 0){
#				break
#			}
	
			activeNW[[2+count]] = activeNW_temp
			activeGVec[allChildren] = 1
			
#			print(activeNW_temp)

			# the active/inactive state of a root node doesn't change with net evolution, as it has no incoming edges
			for (i in nonRootNodes){	
#				print(paste("t ", bla))
#				print(paste("i ", i))
#				print(sum(activeNW_temp[,i]))
#				print(activeNW_temp[,i])
				if ( sum(activeNW_temp[,i]) > 0 & act_mat[i,k]==1) activeGVec[i] = 1
				else activeGVec[i] = 0						
			}
			for (i in nonRootNodes){
				if (activeGVec[i] == 0) activeNW_temp[i,]=0
			}

			geneState[[2+count]][,k] = activeGVec

			parentNodes = which(activeGVec > 0)
			count = count + 1
			
			edgesAdded = unique(edgesAdded)

		} # end while

	} # end of k


	geneStateVec = array(NA,c(n,K,T_))
	T_nw = vector()
	
#	print(activeNW)
	
	for (t in 1:T_){
		colnames(activeNW[[t]]) = rownames(activeNW[[t]]) = nodeNames
		geneStateVec[,,t] = geneState[[t]]
		T_nw = rbind(T_nw, activeNW[[t]])
	}

	T_nw = as.matrix(T_nw)
	
#	print(geneStateVec)
	
	return(list(T_nw=T_nw, geneStateVec=geneStateVec, T_=T_))
}


#source("/home/dx//Thesis/lpNet_Code/src/calcActivation.R")
##source("/home/dx/Dropbox/DX/Thesis/otherCodes/getDdepnNets_inhib.r")

#### data/graph7 from n=10 networks
##n=10
##K=1

##b = c(1,1,1,1,1,1,1,1,1,1,
##			0,1,1,1,1,1,1,1,1,1,
##			1,0,1,1,1,1,1,1,1,1,
##			1,1,0,1,1,1,1,1,1,1,
##			1,1,1,0,1,1,1,1,1,1,
##			1,1,1,1,0,1,1,1,1,1,
##			1,1,1,1,1,0,1,1,1,1,
##			1,1,1,1,1,1,0,1,1,1,
##			1,1,1,1,1,1,1,0,1,1,
##			1,1,1,1,1,1,1,1,0,1,
##			1,1,1,1,1,1,1,1,1,0)

###T_undNet = matrix(c(0,1,0,0,0,0,0,0,0,
###										0,0,-1,0,0,0,0,0,0,
###										0,0,0,0,0,0,0,0,0,
###										0,0,1,0,0,0,0,0,0,
###										0,0,0,1,0,1,0,0,0,
###										0,0,0,0,0,0,1,0,0,
###										0,-1,0,0,0,0,0,1,0,
###										0,0,0,0,0,0,0,0,1,
###										0,0,0,0,0,0,0,0,0), nrow=n,ncol=n,byrow=T)



###T_undNet = matrix(c(0,1,0,0,0,0,0,0,0,
###										0,0,1,0,0,0,0,0,0,
###										0,0,0,1,1,0,0,0,0,
###										0,-1,0,0,0,0,0,0,0,
###										0,0,0,0,0,1,0,0,0,
###										0,0,0,0,0,0,1,0,0,
###										0,0,0,0,0,0,0,1,0,
###										0,0,0,0,0,0,0,0,1,
###										0,0,0,0,0,0,0,0,0), nrow=n,ncol=n,byrow=T)


										
##act_mat = matrix(b, nrow=n, ncol=K)

##act_mat = calcActivation_dyn(T_undNet,b,n,K)									


#n=10
#K=1
###T_undNet = getDDEPNnets_inhib()[[10]]
#bvec = rep(1,n)

##T_undNet = matrix(c(0,1,0,1,0,0,0,0,0,0,
##										1,0,0,0,0,0,1,0,0,0,
##										0,1,0,1,0,0,0,0,1,0,
##										1,0,0,0,0,0,0,-1,0,0,
##										0,0,0,0,0,0,0,1,0,0,
##										0,1,0,1,0,0,0,1,0,0,
##										-1,-1,-1,1,0,0,0,0,0,1,
##										0,0,-1,0,0,0,1,0,0,0,
##										0,0,0,0,0,0,-1,0,0,0,
##										1,0,1,0,0,0,0,0,0,0), nrow=n, ncol=n, byrow=T)


#T_undNet = matrix(c(0,0,0,1,0,0,0,0,0,0,
#										0,0,0,1,0,0,0,1,0,0,
#										0,1,0,0,0,0,0,0,0,0,
#										0,0,0,0,0,0,0,0,1,1,
#										0,1,0,1,0,0,0,0,0,0,
#										0,1,0,0,0,0,0,0,0,0,
#										0,1,0,0,0,0,0,0,0,0,
#										0,0,0,0,0,0,0,0,0,0,
#										1,0,0,0,0,0,1,0,0,0,
#										0,1,0,0,0,0,0,0,0,0), nrow=n,ncol=n,byrow=T)
																				
###bvec = c(1,1,1,1,1,1,1,1,1,1,
###					0,1,1,1,1,1,1,1,1,1,
###					1,0,1,1,1,1,1,1,1,1,
###					1,1,0,1,1,1,1,1,1,1,
###					1,1,1,0,1,1,1,1,1,1,
###					1,1,1,1,0,1,1,1,1,1,
###					1,1,1,1,1,0,1,1,1,1,
###					1,1,1,1,1,1,0,1,1,1,
###					1,1,1,1,1,1,1,0,1,1,
###					1,1,1,1,1,1,1,1,0,1,
###					1,1,1,1,1,1,1,1,1,0)

####print(length(bvec))


#nodeNames = seq(1,n)
#timeSeriesData = build_timeSeriesData(T_undNet, nodeNames,bvec,n,K)


##print("pos chamada")
###T_ = timeSeriesData2$T_
###print(paste("time points ", T_, sep=""))
###print(act_mat)

####save(timeSeriesData, file="net03_geneState.Rdata")
###load(file="net03_geneState.Rdata")

###print(any(timeSeriesData2$geneStateVec == timeSeriesData$geneStateVec) == FALSE)
##print(timeSeriesData$geneStateVec)
###print(timeSeriesData2$geneStateVec)


