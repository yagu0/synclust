#draw (France or...) map with all sites of colors 'cols'
drawMapWithSites = function(M, cols=rep(1,nrow(M)))
{
	xMin = range(M[,1])[1]
	xMax = range(M[,1])[2]
	yMin = range(M[,2])[1]
	yMax = range(M[,2])[2]
	par(mar=c(2,2,2,2))
	plot(0,0,xlim=c(xMin,xMax),ylim=c(yMin,yMax),col="white")
	#plot by color groups (limited to integers)
	maxColor = max(cols)
	for (i in 1:maxColor)
	{
		indices = (1:nrow(M))[cols==i]
		if (length(indices) > 0)
			points(M[indices,1],M[indices,2],col=i,xtitle=NULL)
	}
}

#draw neighborhoods graph on top of a country map (or any other map)
drawNeighborhoodGraph = function(M, NI)
{
	for (i in 1:length(NI))
	{
		for (j in NI[[i]])
			lines(c(M[i,1],M[j,1]),c(M[i,2],M[j,2]))
	}
}

#plot a matrix of curves (in rows)
plotCurves = function(M, cols=rep(1,nrow(M)))
{
	n = nrow(M)
	rg = c(min(M),max(M)) #range for plotting
	for (i in 1:n)
	{
		plot(M[i,],col=cols[i],ylim=rg,type="l")
		if (i < n) par(new=TRUE)
	}
}
