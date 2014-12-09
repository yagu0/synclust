#various util functions for the main program

#preliminary: replace NA's by averaging over each serie's values
#TODO: find a better way to handle missing values
replaceNAs = function(M)
{
	n = nrow(M)
	m = ncol(M)
	res = M
	for (i in 1:n)
	{
		avg = mean(M[i,1:(m-2)] [!is.na(M[i,1:(m-2)])])
		res[i,1:(m-2)] [is.na(M[i,1:(m-2)])] = avg
	}
	return (res)
}

#standardize matrix M (remove mean, divide by standard deviation)
standardize = function(M)
{
	avgM = colMeans(M, na.rm = TRUE)
	stdevs = sqrt( unlist( apply(M, 2, var, na.rm=TRUE) ) )
	res = t(M) - avgM
	res = t(res / stdevs)
	return (list("M"=res,"avg"=avgM,"stv"=stdevs))
}

#opposite of the previous function: get back M from standardized form
destandardize = function(std)
{
	M = std$M
	M = t(M) * std$stv
	M = t(M + std$avg)
	return (M)
}

#remap neighbors into some connex component
remapNeighbors = function(NI, indices)
{
	revIndices = rep(NA, length(NI))
	nc = length(indices)
	for (ii in 1:nc)
		revIndices[ indices[ii] ] = ii
	locNI = list("ix"=as.list(rep(NA,nc)), "ds"=as.list(rep(NA,nc)))
	for (ii in 1:nc)
	{
		locNI$ix[[ii]] = revIndices[ NI$ix[[ indices[ii] ]] ]
		locNI$ds[[ii]] = NI$ds[[ indices[ii] ]]
	}
	return (locNI)
}

#check graph connexity
getConnectedComponents = function(NIix)
{
	return (.Call("getConnectedComponents", NIix));
}

#auxiliary function to display clustering information
promptForMapDisplay = function(stage, coordsM, NIix=NULL, clusters=NULL)
{
	if (is.null(clusters))
		clusters = rep(1, nrow(coordsM))

	shouldDisplay = ""
	if (stage == "interm")
		shouldDisplay = readline(">>> show intermediate map of neighborhoods ? (y/n)\n")
	else if (stage == "final")
	{
		shouldDisplay = readline(
">>> show final map of clusters ? (y/n) \
NOTE: can be plotted later, see '? drawMapWithSites'\n")
	}

	if (shouldDisplay == "y")
	{
		drawMapWithSites(coordsM, clusters)
		if (!is.null(NIix))
			drawNeighborhoodGraph(coordsM,NIix)
		print("Please press 'enter' to continue")
		readline()
		if (!is.null(dev.list()))
			dev.off()
	}
}
