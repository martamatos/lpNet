#
# script to pre-process the experimental datasets (CSV version)
#

options(error = quote({dump.frames(to.file=TRUE); q()}))

args <- commandArgs(TRUE)
dataDir = args[1]

# directory where CSV files are located
setwd(dataDir)


for (CL_number in 1:4){

	sourceNodes = c("c-Met_pY1235", "EGFR_pY1068", "EGFR_pY1173", "EGFR_pY992", "HER2_pY1248", "HER3_pY1298", "FGFR1", "FGFR3")
	cellLines = c("BT20", "BT549", "MCF7", "UACC812")
	inhibitedProteins = list(GSK690693="AKT", GSK690693_GSK1120212="AKT+MEK", PD173074="FGFR1+FGFR3", DMSO="control")

	file01= read.table(file=sprintf("%s_main.csv", cellLines[CL_number]),sep = "\t",header=F)

	file01 = as.matrix(file01)

	inhibitors = unique(file01[,2])
	stimuli = unique(file01[,3])


	if (CL_number != 4){
		antibody_name = unique(file01[2,])
	}

	if (CL_number == 4){
		antibody_name = unique(file01[3,]) 
	}


	timePoints = unique(file01[,4])

	inhibitors = inhibitors[-1]
	inhibitors = inhibitors[-1]
	stimuli = stimuli[-2]
	antibody_name = antibody_name[-1]
	antibody_name = antibody_name[-1]
	timePoints = timePoints[-1]
	timePoints = timePoints[-1]
	timePoints = timePoints[-1]
	timePoints = timePoints[-1]
	if (CL_number == 4) timePoints = timePoints[-1] 

	nInhib = length(inhibitors)
	nStim = length(stimuli)
	nTimePoints = length(timePoints)
	nAntibodies = length(antibody_name)

#	print(inhibitors)
#	print(stimuli)
#	print(antibody_name)
#	print(paste("nAtin", nAntibodies))
#	print(timePoints)
#	print(inhibitedProteins)

	allStim = list()


	nCols = dim(file01)[2]


	for (stim in 1:nStim)
	{
		allStim[[stim]] = array(NA, c(nInhib, nCols, nTimePoints)) 
		rows = which(file01[,3] == stimuli[stim]) # get row numbers corresponding to stim
		tempStim = file01[rows,] # get rows (content) corresponding to stim
		
		for (t in 1:nTimePoints)
		{
			timeRows = which(tempStim[,4] == timePoints[t]) # get row numbers of stim corresponding to time point t
			tempStim_time = tempStim[timeRows,] # get rows (content) of stim corresponding to time point t

			for (inhib in 1:nInhib){

				indRow = which(inhibitors[inhib]==tempStim_time[,2]) # get row numbers of stim + t corresponding to inhib

				if (length(indRow)>0) # if there are no rows, they are not added to the list
				{
					tempInhib=tempStim_time[indRow,]
					if (length(indRow) >= 2) # if there are more than one row per stim/t/inhib average the values (this occurs a t=0min for instance)
					{ 
						mode(tempInhib) = "numeric"
						mean_tempInhib = apply(tempInhib,2,mean, na.rm=T)
						allStim[[stim]][inhib,,t] = c(tempStim_time[indRow,1:4][1,], mean_tempInhib[5:nCols])   
					}
					else allStim[[stim]][inhib,,t] = tempInhib
				}
			} 
		}
	}

	# the first time point is the same for all stimuli
	for (stim in 2:nStim){
		allStim[[stim]][,,1] = allStim[[1]][,,1]
	}
	#print(allStim[[2]])

	#remove low quality proteins

	antibodies_to_remove = c("TAZ_pS89", "FOXO3a_pS318_S321")
	antibodies_to_remove_ind01 = which(antibody_name == antibodies_to_remove[1])
	antibodies_to_remove_ind02 = which(antibody_name == antibodies_to_remove[2])

	if (length(antibodies_to_remove_ind01)>0) antibody_name = antibody_name[-antibodies_to_remove_ind01]
	if (length(antibodies_to_remove_ind02)>0) antibody_name = antibody_name[-antibodies_to_remove_ind02]

	for (stim in 2:nStim){
		if (length(antibodies_to_remove_ind01)>0) allStim[[stim]] = allStim[[stim]][,-antibodies_to_remove_ind01,]
		if (length(antibodies_to_remove_ind02)>0) allStim[[stim]] = allStim[[stim]][,-antibodies_to_remove_ind02,]
	}

	nRows = dim(allStim[[stim]])[2]-4

	# create obs matrix, bvec and delta

	allStimT = list()
	for (stim in 2:nStim)
	{
		nRows = dim(allStim[[stim]])[2]-4 # first 4 columns specify cell line/inhibitor/stimuli/time point
		nCols = dim(allStim[[stim]])[1]
		allStimT[[stim]] = array(NA, c(nRows,nCols, nTimePoints))
		
		# build obs matrix
		mode(allStim[[stim]]) = "numeric"
		for (t in 1:nTimePoints){
			allStimT[[stim]][,,t] = t(allStim[[stim]][,-1:-4,t])
			colnames(allStimT[[stim]]) = inhibitors
			rownames(allStimT[[stim]]) = antibody_name
		}
		
		obs = allStimT[[stim]]	
		
		# normalize data
		for (i in 1:nRows){
			maxI = boxplot(as.vector(obs[i,,]), ylim=c(-3,3), na.rm=T)$stats[5]
			obs[i,,] = obs[i,,]/maxI
		}
		
		
		# build delta and mu
		delta = obs[,which(colnames(obs[,,1])=="DMSO"),1]
		active_mu = rep(NA, nRows)
		active_sd = rep(0.01, nRows)
		inactive_mu = rep(NA, nRows)
		inactive_sd = rep(0.01, nRows)

		for (gene in 1:nRows)
		{
			entries_act = which(obs[gene,,] >= delta[gene], arr.ind=T)
			entries_inact = which(obs[gene,,] < delta[gene], arr.ind=T)
			
			stdGene = sd(as.vector(obs[gene,,]), na.rm=T)
			
			if ( (dim(entries_act)[1] + dim(entries_inact)[1] + length(which(is.na(obs[gene,,])))) != (nCols*nTimePoints)){
			 stop(sprintf("entries act/inact ERROR, expected: %s, actual: %s", nCols*nTimePoints, length(which(is.na(obs[gene, ])))+dim(entries_act)[1]+dim(entries_inact)[1]))
			}
			
			active_ent_mean = c()
			if (dim(entries_act)[1] > 0) {
				for (ent_ac in 1:dim(entries_act)[1]){
					active_ent_mean = c(active_ent_mean, obs[gene, entries_act[ent_ac,1], entries_act[ent_ac,2]])
				}
			}
			else{
				active_ent_mean = c(active_ent_mean, delta[gene]+stdGene)
			}
			
			inactive_ent_mean = c()
			if (dim(entries_inact)[1] > 0) {
				for (ent_inac in 1:dim(entries_inact)[1]){
					inactive_ent_mean = c(inactive_ent_mean, obs[gene, entries_inact[ent_inac,1], entries_inact[ent_inac,2]])
				}
			}
			else{
				inactive_ent_mean = c(inactive_ent_mean, delta[gene]-stdGene)
			}
			
			# safeguard
			if ( !is.na(any(active_ent_mean)) & ((active_ent_mean < delta[gene]) | any (inactive_ent_mean >= delta[gene]))) stop("mu ERROR")
			
			active_mu[gene] = mean(active_ent_mean)
			inactive_mu[gene] = mean(inactive_ent_mean)
		}
		
		# build bvec
		bvec = rep(1, nCols*nRows)
		
		for (inhib in 1:nInhib)
		{	
	
			# get name of the inhibited antibody
			inhibited_proteins_ind = which(names(inhibitedProteins)==inhibitors[inhib])
			inhibited_antibodies = strsplit(inhibitedProteins[[inhibited_proteins_ind]], "\\+")
			
			# get index of antibody, and set respective bvec entry to 0
			for ( inAnt in 1:length(inhibited_antibodies[[1]])){
				inhibited_antibodies_ind=grep(inhibited_antibodies[[1]][inAnt], rownames(obs[,,1]), ignore.case = FALSE, perl=TRUE)
				bvec[(inhib-1)*nRows + inhibited_antibodies_ind] = 0
			}
		}
		
		
		# get the source nodes index
		sourceNodeNumbers = c()
		for ( sNode in 1:length(sourceNodes)){
			sourceNodeNumbers = c(sourceNodeNumbers, grep(sourceNodes[sNode], rownames(obs[,,1]), ignore.case = FALSE, perl=TRUE))
		}

		save(obs, file=sprintf("%s_%s.Rdata", cellLines[CL_number], stimuli[stim]))
		save(delta, file=sprintf("delta_%s_%s.Rdata", cellLines[CL_number], stimuli[stim]))
		save(active_mu, file=sprintf("active_mu_%s_%s.Rdata", cellLines[CL_number], stimuli[stim]))
		save(active_sd, file=sprintf("active_sd_%s_%s.Rdata", cellLines[CL_number], stimuli[stim]))
		save(inactive_mu, file=sprintf("inactive_mu_%s_%s.Rdata", cellLines[CL_number], stimuli[stim]))
		save(inactive_sd, file=sprintf("inactive_sd_%s_%s.Rdata", cellLines[CL_number], stimuli[stim]))
		save(bvec, file=sprintf("bvec_%s_%s.Rdata", cellLines[CL_number], stimuli[stim]))
		save(sourceNodeNumbers, file=sprintf("sourceNodes_%s_%s.Rdata", cellLines[CL_number], stimuli[stim]))
	}
}
