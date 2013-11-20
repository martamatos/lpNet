convertNets_exp = function(cellLines, nCellLines, stimuli, nStim, nAntibodes, function_){

	run=1

	for (cellLine in 1:nCellLines){
		for (stim in 1:nStim){
			setwd(sprintf("01_Code_%s_%s", cellLines[cellLine], stimuli[stim]))
			
			load(sprintf("%s_%s.Rdata", cellLines[cellLine], stimuli[stim]), .GlobalEnv)
			antibodies = rownames(obs)
			annot_node = seq(1,nAntibodes[cellLine])
			
			
			setwd("output")
			dirList=list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
			g=grep(sprintf("%s_run0%s_proc00$", function_, run), dirList, ignore.case = FALSE, perl=TRUE)
			setwd(sprintf("%s", dirList[g]))
			
			load("edges_ddepn.Rdata", .GlobalEnv)

			net = getSampleAdjaMAD_dream8(edges_all=mat,numnodes=nAntibodes[cellLine],annot_node=annot_node,method=median,method2=mad,method2Times=1,septype="->")
			
			rownames(net) = colnames(net) = antibodies
			
			net = net[-dim(net)[1]:-(dim(net)[1]-1), -dim(net)[2]:-(dim(net)[2]-1)]

			setwd("..")
			setwd("..")
			setwd("..")
			
			#convert to sif and eda
			edgesSIF = vector()
			edgesEDA = vector()
			for (i in 1:(nAntibodes[cellLine]-2))
			{
				for (j in 1:(nAntibodes[cellLine]-2)){
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
			
			write(edgesSIF, file=sprintf("lpNet-%s-%s-Network.sif", cellLines[cellLine], stimuli[stim]))
			write(edgesEDA, file=sprintf("lpNet-%s-%s-Network.eda", cellLines[cellLine], stimuli[stim]))
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

cellLines = c("BT20", "BT549", "MCF7", "UACC812")
stimuli = c("EGF", "FGF1", "HGF", "IGF1", "Insulin", "NRG1", "PBS", "Serum")
nAntibodes = c(48, 46, 41, 46)
nCellLines = length(cellLines)
nStim = length(stimuli)

convertNets_exp(cellLines, nCellLines, stimuli, nStim, nAntibodes, function_)

