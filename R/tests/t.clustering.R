#test several clustering methods on iris dataset (setosa should be found)
test.clustering1 = function()
{
	data(iris)

	#get neighborhoods from data [25 is high, but shouldn't be lower to have 1 connex comp.]
	NI = .Call("getNeighbors", as.matrix(iris[,1:4]), 25, 0.0, 1, TRUE)

	for (dtype in c("spath"))#,"ectd")) #bug: TODO
	{
		#get distances from neighborhoods; should be OK for all distances 
		#except "simple" (which is treated as a special case with built-in R funcs)
		distances = synclust:::getDistances(dtype, NI)

		for (cmeth in c("KM","HC"))
		{
			#finally, get clusters
			clusters = synclust:::getClusters(distances, cmeth, K=3)
			#check that cluster 'setosa' is pure and separated
			uqclust = unique(clusters[1:50])
			checkTrue(length(uqclust) == 1)
			checkTrue(! uqclust[1] %in% clusters[51:150])
		}
	}
}

#test several parameters agencements on custom non-isotropic gaussian dataset (2D)
test.clustering2 = function()
{
	clustSize = 33
	
	require(mvtnorm)
	set.seed(32)
	gaussian1 = rmvnorm(clustSize, mean=c(-4.0,-6.0), sigma=matrix(c(1.0,0.7,0.7,1.0),nrow=2))
	gaussian2 = rmvnorm(clustSize, mean=c(0.0,0.0), sigma=matrix(c(1.0,0.0,0.0,1.0),nrow=2))
	gaussian3 = rmvnorm(clustSize, mean=c(4.0,-6.0), sigma=matrix(c(1.0,-0.7,-0.7,1.0),nrow=2))
	M = rbind(gaussian1, rbind(gaussian2, gaussian3))
	
	#get neighborhoods from data [25 is high, but shouldn't be much lower to have 1 connex comp.]
	NI = .Call("getNeighbors", M, 25, 0.0, 1, TRUE)
	
	for (dtype in c("spath"))#,"ectd")) #TODO
	{
		#get distances from neighborhoods; should be OK for all distances 
		#except "simple" (which is treated as a special case with built-in R funcs)
		distances = synclust:::getDistances(dtype, NI)
		
		for (cmeth in c("KM","HC"))
		{
			#finally, get clusters
			clusters = synclust:::getClusters(distances, cmeth, K=3)
			
			#soft check, because have to succeed at each run
			srt = sort(clusters)
			checkTrue( sum( srt[1:clustSize] == 1 ) / clustSize >= 0.8 )
			checkTrue( sum( srt[(clustSize+1):(2*clustSize)] == 2 ) / clustSize >= 0.8 )
			checkTrue( sum( srt[(2*clustSize+1):(3*clustSize)] == 3 ) / clustSize >= 0.8 )
		}
	}
}

#test several parameters agencements on custom "two moons one circle" dataset (2D)
test.clustering3 = function()
{
	clustSize = 150

	set.seed(32)	
	M = matrix(nrow=3*clustSize,ncol=2)
	#big circle: radius = 10
	rdata = runif(clustSize, min=0, max=2*pi)
	M[1:clustSize,] = 10 * cbind(cos(rdata), sin(rdata))
	#moon 1: half circle of radius 5 centered at (-2, -0.5)
	rdata = runif(clustSize, min=0, max=pi)
	M[(clustSize+1):(2*clustSize),] = cbind(5*cos(rdata)-2, 5*sin(rdata)-0.5)
	#moon 2: half circle of radius 5 centered at (2, 0.5)
	rdata = runif(clustSize, min=pi, max=2*pi)
	M[(2*clustSize+1):(3*clustSize),] = cbind(5*cos(rdata)+2, 5*sin(rdata)+0.5)
	
	#add small global noise
	M = M + rnorm(2*clustSize,sd=0.1)
	
	#get neighborhoods from data [25 is high, but shouldn't be much lower to have 1 connex comp.]
	NI = .Call("getNeighbors", M, 25, 0.0, 1, TRUE)
	
	#only ECTD distance can be used, because forcing connexity implies 
	#creating shortcuts in graph, which strongly affect spath distance
	distances = synclust:::getDistances("ectd", NI)
	
	#only hierarchical clustering can succeed here
	clusters = synclust:::getClusters(distances, "HC", K=3)

	srt = sort(clusters)
	checkTrue( sum( srt[1:clustSize] == 1 ) / clustSize >= 0.90 )
	checkTrue( sum( srt[(clustSize+1):(2*clustSize)] == 2 ) / clustSize >= 0.90 )
	checkTrue( sum( srt[(2*clustSize+1):(3*clustSize)] == 3 ) / clustSize >= 0.90 )
}

#renumbering if clusters have too high labels
test.reordering = function()
{
	clusters = c(1,6,8,8,8,1,1,1,6,6,6,8,8,1,1,6,8)
	checkEquals(sort(unique(synclust:::reordering(clusters))),c(1,2,3))
	clusters = c(23,3,23,77,77,77,1,12,12,12,77,12,23,23,12,23,77,12,23,77,1)
	checkEquals(sort(unique(synclust:::reordering(clusters))),c(1,2,3,4,5))
}
