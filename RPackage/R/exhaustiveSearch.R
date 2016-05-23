exhaustiveSearch <- function(graph, interestVertices, countDisconnected)
{
	if(class(graph) == "igraph")
	{
		if(igraph::is.directed(graph))
		{
			stop("Input `graph' must be undirected")
		}
		start <- Sys.time()
		result <- .Call("exhaustiveSearch_igraph", graph, interestVertices, countDisconnected, PACKAGE="networkReliability")
		end <- Sys.time()
	}
	else if(class(graph) == "graphNEL")
	{
		start <- Sys.time()
		result <- .Call("exhaustiveSearch_graphNEL", graph, interestVertices, countDisconnected, PACKAGE="networkReliability")
		end <- Sys.time()
	}
	else if(class(graph) == "graphAM")
	{
		start <- Sys.time()
		result <- .Call("exhaustiveSearch_graphAM", graph, interestVertices, countDisconnected, PACKAGE="networkReliability")
		end <- Sys.time()
	}
	else
	{
		stop("Input graph must have class \"igraph\", \"graphAM\" or \"graphNEL\"")
	}
	return(new("exhaustiveSearchResult", data = mpfr(result), call = match.call(), start = start, end = end, interestVertices = as.integer(interestVertices), graph = graph, countDisconnected = countDisconnected))

}
