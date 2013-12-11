ddepn_realTime = function(obs_modified)
{
	obs_modified = round(obs_modified, digits=5)
	
	delta = obs_modified[,,1]
	delta[3:14,]=delta[3:14,]+0.00001
	delta = round(delta, digits=5)

	active_mu = matrix(NA,nrow=16, ncol=3)
	active_sd = matrix(0.01,nrow=16, ncol=3)
	inactive_mu = matrix(NA,nrow=16, ncol=3)
	inactive_sd = matrix(0.01,nrow=16, ncol=3)

	for (e in 1:K){
		for (i in 1:n){ 
			vec = obs_modified[i,e,]
			active_mu[i,e] =  mean(vec[which(vec >= delta[i,e])], na.rm=T)
			inactive_mu[i,e] = mean(vec[which(vec < delta[i,e])], na.rm=T)
			
		}
	}
	
	for (i in 1:n){
		active_mu[i,which(is.na(active_mu[i,]))]=mean(active_mu[i,], na.rm=T)
		inactive_mu[i,which(is.na(inactive_mu[i,]))]=mean(inactive_mu[i,], na.rm=T)
	}
	
	return(list(obs_modified=obs_modified, delta=delta, active_mu=active_mu, active_sd=active_sd, inactive_mu=inactive_mu, inactive_sd=inactive_sd))
}


dream_insilico_realTime = function(obs_modified)
{
	nStim = 8
	nInhib = 4
	nAntibodies = 20
	sourceNodes = c("AB10", "AB12", "AB20")
	sNodes_num = c(10,12,20)
	nRows = dim(obs_modified)[1]
	nTimePoints = dim(obs_modified)[3]
	obs_modified = round(obs_modified, digits=2)
	
	delta_temp = obs_modified[,which(colnames(obs_modified[,,1])=="None"),1]
	delta = vector()
	print(delta_temp)
	print(nStim)
	for (stim in 1:nStim){
		delta = cbind(delta, matrix(rep(delta_temp[,stim], 4), nrow=nAntibodies, ncol=nInhib))
	}
	stop("A")
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
			entries_act = which(obs_modified[gene,((stim-1)*nInhib+1):((stim-1)*nInhib+nInhib),] >= delta[gene, ((stim-1)*nInhib+1)])
			entries_inact = which(obs_modified[gene,((stim-1)*nInhib+1):((stim-1)*nInhib+nInhib),] < delta[gene, ((stim-1)*nInhib+1)])
			
			
			if ( (length(entries_act)+ length(entries_inact) + length(which(is.na(obs_modified[gene,((stim-1)*nInhib+1):((stim-1)*nInhib+nInhib),])))) != (nInhib*nTimePoints)){
			 stop(sprintf("entries act/inact ERROR, expected: %s, actual: %s", nInhib*nTimePoints, length(which(is.na(obs_modified[gene,((stim-1)*nInhib+1):((stim-1)*nInhib+nInhib),]))) + length(entries_act)+length(entries_inact)))
			}
			
			stdGene = sd(as.vector(obs_modified[gene,((stim-1)*nInhib+1):((stim-1)*nInhib+nInhib),]), na.rm=T)

			active_ent_mean = c()
			if (length(entries_act) > 0) {
				active_ent_mean = obs_modified[gene,((stim-1)*nInhib+1):((stim-1)*nInhib+nInhib),][entries_act]
			}

			inactive_ent_mean= c()
			if (length(entries_inact) > 0) {
				inactive_ent_mean = obs_modified[gene,((stim-1)*nInhib+1):((stim-1)*nInhib+nInhib),][entries_inact]
			}

			# safeguard
			if (any(active_ent_mean < delta[gene, ((stim-1)*nInhib+1)]) | any (inactive_ent_mean >= delta[gene, ((stim-1)*nInhib+1)])) stop("mu ERROR")
			
			active_mu[gene, ((stim-1)*nInhib+1):((stim-1)*nInhib+nInhib)] = mean(active_ent_mean)
			inactive_mu[gene, ((stim-1)*nInhib+1):((stim-1)*nInhib+nInhib)] = mean(inactive_ent_mean)
		}
	}

	active_mu = round(active_mu, digits=2)
	inactive_mu = round(inactive_mu, digits=2)
	
	return(list(obs_modified=obs_modified, delta=delta, active_mu=active_mu, active_sd=active_sd, inactive_mu=inactive_mu, inactive_sd=inactive_sd))
}
