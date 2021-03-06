approximateZeroVarianceWOR <- function(graph, probability, n, seed, interestVertices)
{
	if(class(graph) == "igraph")
	{
		if(igraph::is.directed(graph))
		{
			stop("Input `graph' must be undirected")
		}
		start <- Sys.time()
		result <- .Call("approximateZeroVarianceWOR_igraph", graph, probability, n, seed, interestVertices, PACKAGE="networkReliability")
		end <- Sys.time()
	}
	else if(class(graph) == "graphNEL")
	{
		start <- Sys.time()
		result <- .Call("approximateZeroVarianceWOR_graphNEL", graph, probability, n, seed, interestVertices, PACKAGE="networkReliability")
		end <- Sys.time()
	}
	else if(class(graph) == "graphAM")
	{
		start <- Sys.time()
		result <- .Call("approximateZeroVarianceWOR_graphAM", graph, probability, n, seed, interestVertices, PACKAGE="networkReliability")
		end <- Sys.time()
	}
	else
	{
		stop("Input graph must have class \"igraph\", \"graphAM\" or \"graphNEL\"")
	}
	return(new("approximateZeroVarianceWORResult", estimate = mpfr(result), call = match.call(), start = start, end = end, n = as.integer(n), interestVertices = as.integer(interestVertices), seed = as.integer(seed), graph = graph, probability = probability))
}
