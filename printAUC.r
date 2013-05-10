printAUC = function(auc,sd_all,randROC,totalruns, t, ouputDir)
{	

#		print("randROC")
#		print(randROC)
		auc2 = vector()
		auc2 = rbind(auc2, auc)
		auc2 = rbind(auc2, c(mean(randROC),sd(randROC),NA,NA))
		
#		print("auc2 mean")
#		print(auc2)
		
#		auc2 = rbind(auc2, c(randROC[[2]],NA,NA,NA))

		colnames(auc2) <- c("mean AUC-ROC","mean AUC-PR","sd AUC-ROC","sd AUC-PR")
		rownames(auc2) = c(sd_all, "AUC-ROC_rand")
		
		write.table(auc2,sprintf('%s/auc_t=%s.txt', ouputDir,t),quote=FALSE,row.names=F)
	

#		postscript(sprintf('%s/aucROC_t=%s.eps', ouputDir,t),horizontal=T)
#		plot(sd_all,auc[1:3,1],main="AUC values for varying noise",xlab="Noise",ylab="AUC value",ylim=c(0,1),type="l",lwd=2)
		
#		lines(sd_all,auc[,1]+(auc[,3])/sqrt(totalruns),col="black",lty="dotted",lwd=1)
#		lines(sd_all,auc[,1]-(auc[,3])/sqrt(totalruns),col="black",lty="dotted",lwd=1)
		
#		lines(sd_all,auc[,2],col="red",lty=1,lwd=2)
		
#		lines(sd_all,auc[,2]+(auc[,4])/sqrt(totalruns),col="red",lty="dotted",lwd=1)
#		lines(sd_all,auc[,2]-(auc[,4])/sqrt(totalruns),col="red",lty="dotted",lwd=1)
		
#		abline(h=randROC[1],col="black",lty=3,lwd=1.5)
#		abline(h=randROC[2],col="red",lty=3,lwd=1.5)
		
#		legend("topright",legend=c("AU-ROC curve","AU-PR curve","AU-ROC random","PR-ROC random"),col=c("black","red","black","red"),lty=c(1,1,3,3),lwd=c(1.5,1.5,1.5,1.5))
		
#		dev.off()
}


printAUC_und = function(auc,sd_all,randROC,randPR,totalruns, t, ouputDir)
{
#		print("randROC")
#		print(randROC)
		auc2 = sd_all
		auc2 =  cbind(auc2, auc)
		
#		auc_rand = vector()
#		temp = c(mean(randROC),mean(randPR),sd(randROC),sd(randPR))
		
#		for (it in 1:length(sd_all)){
#			auc_rand = rbind(auc_rand,c(sd_all[it], temp))
#		}
		
		
		
#		print("auc2 mean")
#		print(auc2)
		


#		colnames(auc2) <- c("std","mean AUC-ROC","mean AUC-PR","sd AUC-ROC","sd AUC-PR")
#		rownames(auc2) = c(sd_all)
		
##		print(auc2)
#		print(getwd())
		write.table(auc2,sprintf('%s/auc_und.txt', ouputDir),quote=FALSE,row.names=F,col.names=F)
#		write.table(t(auc_rand),sprintf('%s/auc_und_randROC.txt', ouputDir),quote=FALSE,row.names=F,col.names=F)
	

#		postscript(sprintf('%s/aucROC_t=%s.eps', ouputDir,t),horizontal=T)
#		plot(sd_all,auc[1:3,1],main="AUC values for varying noise",xlab="Noise",ylab="AUC value",ylim=c(0,1),type="l",lwd=2)
		
#		lines(sd_all,auc[,1]+(auc[,3])/sqrt(totalruns),col="black",lty="dotted",lwd=1)
#		lines(sd_all,auc[,1]-(auc[,3])/sqrt(totalruns),col="black",lty="dotted",lwd=1)
		
#		lines(sd_all,auc[,2],col="red",lty=1,lwd=2)
		
#		lines(sd_all,auc[,2]+(auc[,4])/sqrt(totalruns),col="red",lty="dotted",lwd=1)
#		lines(sd_all,auc[,2]-(auc[,4])/sqrt(totalruns),col="red",lty="dotted",lwd=1)
		
#		abline(h=randROC[1],col="black",lty=3,lwd=1.5)
#		abline(h=randROC[2],col="red",lty=3,lwd=1.5)
		
#		legend("topright",legend=c("AU-ROC curve","AU-PR curve","AU-ROC random","PR-ROC random"),col=c("black","red","black","red"),lty=c(1,1,3,3),lwd=c(1.5,1.5,1.5,1.5))
		
#		dev.off()
}
