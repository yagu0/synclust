#main function (choice between kmeans and hierarchical clustering)
getClusters = function(distances, method, K)
{
	clusts = c()
	if (method=="KM")
	{
		nstart = 10 #number of kmeans random restarts
		maxiter = 100 #maximum iterations count in each km run
		clusts = .Call("kmeansWithDistances", distances, K, nstart, maxiter)
	}
	else if (method=="HC")
	{
		#simple hierarchical clustering using ECT distances
		hct = hclust(as.dist(distances),method="ward.D")
		clusts = cutree(hct, K)
	}
	return (clusts)
}

# renumbering step (post-processing after clustering)
reordering = function(clusts)
{
	newCl = clusts
	maxInd = max(clusts)
	counter = 1
	for (i in 1:maxInd)
	{
		if (sum(clusts == i) > 0)
		{
			newCl[clusts == i] = counter
			counter = counter + 1
		}
	}
	return (newCl)
}
