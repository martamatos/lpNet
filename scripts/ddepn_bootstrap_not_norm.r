setwd("../data")
load("hcc1954.RData")

#----------------------------------------------------------------------
# treat data
#----------------------------------------------------------------------

#sumData = list()
#ind=0

#for (e in 1:3){
#  sumData[[e]] = vector()
#  for (i in 1:10){  
#    sumData[[e]] = cbind(sumData[[e]], apply(hcc1954[,(ind+1):(ind+15)], 1, mean ) )
#    ind = ind+15
#  }
#}

#print(sumData)
#for (e in 1:3){
#  for (i in 1:16){
#    sumData[[e]][i,] =  sumData[[e]][i,]/boxplot(sumData[[e]][i,], na.rm=T)$stats[5,]
#  }
#}

ind = 0
dataFinal = array(NA,c(16,3,10,15))
for (e in 1:3){
	for (t in 1:10){
		dataFinal[,e,t,] = hcc1954[,(ind+1):(ind+15)]
		ind = ind+15
	}
}




dataFinal = round(dataFinal, digits=1)
print(dataFinal)
#save(dataFinal, file="dataFinal.Rdata")



#----------------------------------------------------------------------
# get delta and in/active_mu/sd
#----------------------------------------------------------------------


delta = apply(dataFinal[,,1,], c(1,2), mean)

print("------------")
print(delta)
delta[3:14,]=delta[3:14,]+0.1
print(delta)
print("------------")

delta = round(delta, digits=1)

active_mu = matrix(NA,nrow=16, ncol=3)
active_sd = matrix(0.01,nrow=16, ncol=3)
inactive_mu = matrix(NA,nrow=16, ncol=3)
inactive_sd = matrix(0.01,nrow=16, ncol=3)

for (exp in 1:3){
  for (i in 1:16){ 

    vec = dataFinal[i,exp,,]
#     print(vec)
    active_mu[i,exp] =  mean(vec[which(vec >= delta[i,exp])], na.rm=T)
    inactive_mu[i,exp] = mean(vec[which(vec < delta[i,exp])], na.rm=T)
    
  }
}


for (i in 1:16){
  active_mu[i,which(is.na(active_mu[i,]))]=mean(active_mu[i,], na.rm=T)
  inactive_mu[i,which(is.na(inactive_mu[i,]))]=mean(inactive_mu[i,], na.rm=T)
}

active_mu = round(active_mu, digits=1)
inactive_mu = round(inactive_mu, digits=1)


# inactive_mu[which(is.na(inactive_mu))]=0
print(delta)
print(active_mu)
print(inactive_mu)

save(delta, file="delta_ddepn14_not_norm.Rdata")
save(active_mu, file="active_mu_ddepn14_not_norm.Rdata")
save(active_sd, file="active_sd_ddepn14_not_norm.Rdata")
save(inactive_mu, file="inactive_mu_ddepn14_not_norm.Rdata")
save(inactive_sd, file="inactive_sd_ddepn14_not_norm.Rdata")

save(dataFinal, file="dataFinal14_not_norm.Rdata")

