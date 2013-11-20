convertNets_insilico = function(cellLines, nCellLines, stimuli, nStim, nAntibodes, function_){

	run=1
	
	for (cellLine in 1:nCellLines){
		for (stim in 1:nStim){
			setwd(sprintf("01_Code_%s", cellLines[cellLine]))
			
			load(sprintf("%s.Rdata", cellLines[cellLine]), .GlobalEnv)
			antibodies = rownames(obs)
			annot_node = seq(1,nAntibodes[cellLine])
			
			setwd("output")
			dirList=list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
			g=grep(sprintf("%s_run0%s_proc00$", function_, run), dirList, ignore.case = FALSE, perl=TRUE)
			setwd(sprintf("%s", dirList[g]))
			
			load("edges_ddepn.Rdata", .GlobalEnv)

			net = getSampleAdjaMAD_dream8(edges_all=mat,numnodes=nAntibodes[cellLine],annot_node=annot_node,method=median,method2=mad,method2Times=1,septype="->")

			rownames(net) = colnames(net) = antibodies
						
			setwd("..")
			setwd("..")
			setwd("..")
			
			#convert to sif and eda
			edgesSIF = vector()
			edgesEDA = vector()
			for (i in 1:nAntibodes[cellLine])
			{
				for (j in 1:nAntibodes[cellLine]){
					if (net[i,j] !=0){
						if (net[i,j] < 0){
							edgesSIF= c(edgesSIF, sprintf("%s -1 %s", rownames(net)[i], colnames(net)[j]))
							edgesEDA = c(edgesEDA, sprintf("%s (-1) %s = %s", rownames(net)[i], colnames(net)[j], abs(net[i,j])))	
						}
						if (net[i,j] > 0){
							edgesSIF= c(edgesSIF, sprintf("%s 1 %s", rownames(net)[i], colnames(net)[j]))
							edgesEDA = c(edgesEDA, sprintf("%s (1) %s = %s", rownames(net)[i], colnames(net)[j], abs(net[i,j])))

						}
					}
				}
			}
			
			write(edgesSIF, file="lpNet-Insilico-Network.sif")
			write(edgesEDA, file="lpNet-Insilico-Network.eda")
		}
	}

}


args <- commandArgs(TRUE)
dataDir = args[1]

source("lpNet-Code/sourceDir.R")
sourceDir("lpNet-Code")


# directory where results are
print(dataDir)
setwd(dataDir)

function_="doILP_dyn_discretized_dream8_new"

cellLines = c("insilico")
stimuli = c("")
nAntibodes = c(20)
nCellLines = length(cellLines)
nStim = length(stimuli)

convertNets_insilico(cellLines, nCellLines, stimuli, nStim, nAntibodes, function_)

