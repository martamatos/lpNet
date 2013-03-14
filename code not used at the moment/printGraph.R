printGraph = function(res, n, i, graphSize_)
{
		adj=getAdja(res,n)
		write.matrix(adj, file = paste("output/adj",i, sep=""), sep = " ")
		#print(paste("output/adj",i, sep=""))
		#print(sprintf('python net.py  printGraph_%s output/adj%s output/graph%s', graphSize_, i, i))
		system(sprintf('python net.py  printGraph_%s output/adj%s output/graph%s', graphSize_, i, i), intern=TRUE)
		print(adj)
}
