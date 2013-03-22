printAUC_median = function(auc,sd_all,randROC,totalruns, t, ouputDir)
{
		auc2 = vector()
		auc2 = rbind(auc2, auc)
		auc2 = rbind(auc2, c(quantile(randROC),NA,NA,NA,NA,NA))
		print("auc2 median")
		print(auc2)
#		auc2 = rbind(auc2, c(randROC[[2]],NA,NA,NA,NA,NA,NA,NA,NA,NA))

		colnames(auc2) <- c("quant1 AUC-ROC","quant2 AUC-ROC","quant3 AUC-ROC","quant4 AUC-ROC","quant5 AUC-ROC",
												"quant1 AUC-PR","quant2 AUC-PR","quant3 AUC-PR","quant4 AUC-PR","quant5 AUC-PR")
		rownames(auc2) = c(sd_all, "AUC-ROC_rand")
		
#		print("ola")
#		print(auc2)
#		print(ouputDir)
#		print(t)
		write.table(auc2,sprintf('%s/auc_median_t=%s.txt', ouputDir,t),quote=FALSE,row.names=F)
	

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



printAUC_median_und = function(auc,sd_all,randROC,randPR,totalruns, t, ouputDir)
{
		auc2 = sd_all
		auc2 = cbind(auc2, auc)
		auc2 = rbind(auc2, c(0,quantile(randROC),quantile(randPR)))
		
		
#		print("auc2 median")
#		print(auc2)
#		auc2 = rbind(auc2, c(randROC[[2]],NA,NA,NA,NA,NA,NA,NA,NA,NA))

#		colnames(auc2) <- c("quant1 AUC-ROC","quant2 AUC-ROC","quant3 AUC-ROC","quant4 AUC-ROC","quant5 AUC-ROC",
#												"quant1 AUC-PR","quant2 AUC-PR","quant3 AUC-PR","quant4 AUC-PR","quant5 AUC-PR")
#		rownames(auc2) = c(sd_all, "AUC-ROC_rand")
		
#		print(auc2)
		write.table(auc2,sprintf('%s/auc_und_median.txt', ouputDir),quote=FALSE,row.names=F,col.names=F)
	

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
