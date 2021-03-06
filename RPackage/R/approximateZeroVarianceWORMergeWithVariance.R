approximateZeroVarianceWORMergeWithVariance <- function(graph, probability, n, seed, interestVertices, graphFile = "")
{
	if(class(graph) == "igraph")
	{
		if(igraph::is.directed(graph))
		{
			stop("Input `graph' must be undirected")
		}
		start <- Sys.time()
		result <- .Call("approximateZeroVarianceWORMergeWithVariance_igraph", graph, probability, n, seed, interestVertices, graphFile, PACKAGE="networkReliability")
		end <- Sys.time()
	}
	else if(class(graph) == "graphNEL")
	{
		start <- Sys.time()
		result <- .Call("approximateZeroVarianceWORMergeWithVariance_graphNEL", graph, probability, n, seed, interestVertices, graphFile, PACKAGE="networkReliability")
		end <- Sys.time()
	}
	else if(class(graph) == "graphAM")
	{
		start <- Sys.time()
		result <- .Call("approximateZeroVarianceWORMergeWithVariance_graphAM", graph, probability, n, seed, interestVertices, graphFile, PACKAGE="networkReliability")
		end <- Sys.time()
	}
	else
	{
		stop("Input graph must have class \"igraph\", \"graphAM\" or \"graphNEL\"")
	}
	return(new("approximateZeroVarianceWORWithVarianceResult", estimate = mpfr(result$estimate), varianceEstimate = mpfr(result$varianceEstimate), call = match.call(), start = start, end = end, n = as.integer(n), interestVertices = as.integer(interestVertices), seed = as.integer(seed), graph = graph, probability = probability))
}
