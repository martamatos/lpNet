#
# script to pre-process the insilico dataset (CSV version)
#

options(error = quote({dump.frames(to.file=TRUE); q()}))

args <- commandArgs(TRUE)
dataDir = args[1]
outputDir = args[2]

# directory where the insilico.csv is located
setwd(dataDir)

CL_number = 1

sourceNodes = c("AB10", "AB12", "AB20")
sNodes_num = c(10,12,20)

cellLines = c("insilico")
inhibitedProteins = list(INH1="AB12", INH2="AB5", INH3="AB8")

file01 = read.table(file="insilico.csv",sep = "\t",header=F)


file01 = as.matrix(file01)

inhibitors = unique(file01[,2])
stimuli = unique(file01[,3])
antibody_name = unique(file01[1,])

timePoints = unique(file01[,4])

#inhibitors = inhibitors[-1]
inhibitors = inhibitors[-1]
stimuli = stimuli[-1]
antibody_name = antibody_name[-1]
antibody_name = antibody_name[-1]
antibody_name = antibody_name[-1]
antibody_name = antibody_name[-1]
timePoints = timePoints[-1]
timePoints = timePoints[-length(timePoints)]

nInhib = length(inhibitors)
nStim = length(stimuli)
nTimePoints = length(timePoints)
nAntibodies = length(antibody_name)

#print(inhibitors)
#print(stimuli)
#print(antibody_name)
#print(nAntibodies)
#print(timePoints)
#print(inhibitedProteins)

allStim = list()


nCols = dim(file01)[2]


for (stim in 1:nStim)
{
#  allStim[[stim]] = array(NA, c(nInhib, nCols, nTimePoints, nReplicates))
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
        print(tempInhib)
        if (length(indRow) >= 2) # if there are more than one row per stim/t/inhib average the values (this occurs a t=0min for instance)
        { 
          mode(tempInhib) = "numeric"
          mean_tempInhib = apply(tempInhib,2,mean)
          allStim[[stim]][inhib,,t] = c(tempStim_time[indRow,1:4][1,], mean_tempInhib[5:nCols])   
        }
        else allStim[[stim]][inhib,,t] = tempInhib
#        print(tempInhib)
#        allStim[[stim]][inhib,,t] = tempInhib
      }
    } 
  }
}


print(allStim)
# create obs matrix, bvec and delta

allStimT = list()
bvec = c()
obs_list = list()
for (t in 1:nTimePoints){
	obs_list[[t]] = vector()
}


for (stim in 1:nStim)
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
	obs_temp = allStimT[[stim]]
	
	# build bvec
	bvec_temp = rep(1, nCols*nRows)
	
	for (inhib in 2:nInhib)
	{	
		# get name of the inhibited antibody
		inhibited_proteins_ind = which(names(inhibitedProteins)==inhibitors[inhib])
		inhibited_antibodies = strsplit(inhibitedProteins[[inhibited_proteins_ind]], "\\+")
				
		# get index of antibody, and set respective bvec entry to 0
		for ( inAnt in 1:length(inhibited_antibodies[[1]])){
			inhibited_antibodies_ind=grep(inhibited_antibodies[[1]][inAnt], rownames(obs_temp[,,1]), ignore.case = FALSE, perl=TRUE)
			bvec_temp[(inhib-1)*nRows + inhibited_antibodies_ind] = 0
		}
	}

	bvec = c(bvec, bvec_temp)
	for (t in 1:nTimePoints){
		obs_list[[t]] = cbind(obs_list[[t]], obs_temp[,,t])
	}
}

obs = array(NA, c(nAntibodies, nInhib*nStim, nTimePoints))
for (t in 1:nTimePoints){
	obs[,,t] = obs_list[[t]]
}

rownames(obs) = antibody_name
colnames(obs) = rep(inhibitors, nStim)


# build delta
delta_temp = obs[,which(colnames(obs[,,1])=="None"),1]
delta = vector()

for (stim in 1:nStim){
	delta = cbind(delta, matrix(rep(delta_temp[,stim], 4), nrow=nAntibodies, ncol=nInhib))
}

obs = round(obs, digits=2)
delta = round(delta, digits=2)

for (i in 1:nAntibodies){
	if (!(i %in% sNodes_num)){
		delta[i,] = delta[i,] + 0.01
	}
}

active_mu = matrix(NA, nrow=nAntibodies, ncol=nStim*nInhib)
active_sd = matrix(0.01, nrow=nAntibodies, ncol=nStim*nInhib)
inactive_mu = matrix(NA, nrow=nAntibodies, ncol=nStim*nInhib)
inactive_sd = matrix(0.01, nrow=nAntibodies, ncol=nStim*nInhib)


for (gene in 1:nRows)
{
	for (stim in 1:nStim){
		entries_act = which(obs[gene,((stim-1)*nInhib+1):((stim-1)*nInhib+nInhib),] >= delta[gene, ((stim-1)*nInhib+1)])
		entries_inact = which(obs[gene,((stim-1)*nInhib+1):((stim-1)*nInhib+nInhib),] < delta[gene, ((stim-1)*nInhib+1)])
		
		
		if ( (length(entries_act)+ length(entries_inact) + length(which(is.na(obs[gene,((stim-1)*nInhib+1):((stim-1)*nInhib+nInhib),])))) != (nInhib*nTimePoints)){
		 stop(sprintf("entries act/inact ERROR, expected: %s, actual: %s", nInhib*nTimePoints, length(which(is.na(obs[gene,((stim-1)*nInhib+1):((stim-1)*nInhib+nInhib),]))) + length(entries_act)+length(entries_inact)))
		}
		
		stdGene = sd(as.vector(obs[gene,((stim-1)*nInhib+1):((stim-1)*nInhib+nInhib),]), na.rm=T)

		
		active_ent_mean = c()
		if (length(entries_act) > 0) {
			active_ent_mean = obs[gene,((stim-1)*nInhib+1):((stim-1)*nInhib+nInhib),][entries_act]
		}

		inactive_ent_mean= c()
		if (length(entries_inact) > 0) {
			inactive_ent_mean = obs[gene,((stim-1)*nInhib+1):((stim-1)*nInhib+nInhib),][entries_inact]
		}

		# safeguard
		if (any(active_ent_mean < delta[gene, ((stim-1)*nInhib+1)]) | any (inactive_ent_mean >= delta[gene, ((stim-1)*nInhib+1)])) stop("mu ERROR")
		
		active_mu[gene, ((stim-1)*nInhib+1):((stim-1)*nInhib+nInhib)] = mean(active_ent_mean)
		inactive_mu[gene, ((stim-1)*nInhib+1):((stim-1)*nInhib+nInhib)] = mean(inactive_ent_mean)
	}
}


sourceNodeNumbers = c()
for ( sNode in 1:length(sourceNodes)){
	sourceNodeNumbers = c(sourceNodeNumbers, grep(sourceNodes[sNode], rownames(obs[,,1]), ignore.case = FALSE, perl=TRUE))
}

active_mu = round(active_mu, digits=2)
inactive_mu = round(inactive_mu, digits=2)


setwd(outputDir)

save(obs, file=sprintf("%s.Rdata", cellLines[CL_number]))
save(delta, file=sprintf("delta_%s.Rdata", cellLines[CL_number]))
save(active_mu, file=sprintf("active_mu_%s.Rdata", cellLines[CL_number]))
save(active_sd, file=sprintf("active_sd_%s.Rdata", cellLines[CL_number]))
save(inactive_mu, file=sprintf("inactive_mu_%s.Rdata", cellLines[CL_number]))
save(inactive_sd, file=sprintf("inactive_sd_%s.Rdata", cellLines[CL_number]))
save(bvec, file=sprintf("bvec_%s.Rdata", cellLines[CL_number]))

