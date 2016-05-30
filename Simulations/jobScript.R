source("./generateScenarios.R")
SCENARIO_INDEX <- as.integer(Sys.getenv("SCENARIO_INDEX"))
cat("SCENARIO_INDEX=", SCENARIO_INDEX, "\n", sep="")
nReps <- scenarios[SCENARIO_INDEX, "nReps"]

outputFile <- file.path("results", scenarios[SCENARIO_INDEX, "file"])
tmpFile <- paste0(outputFile, ".tmp")
library(networkReliability)
library(stringr)
library(igraph)

graph <- scenarios[SCENARIO_INDEX, "graph"]
method <- scenarios[SCENARIO_INDEX, "method"]
sampleSize <- scenarios[SCENARIO_INDEX, "sampleSize"]
graph <- scenarios[SCENARIO_INDEX, "graph"]
probability <- scenarios[SCENARIO_INDEX, "probability"]
nReps <- scenarios[SCENARIO_INDEX, "nReps"]

if(graph == "dodecahedron")
{
	graph <- igraph::read_graph(file = system.file("data/dodecahedron.graphml", package="networkReliability"), format = "graphml")
	interestVertices <- c(1, 20)
} else if(graph == "grid9")
{
	graph <- igraph::make_lattice(dim = 2, length = 9)
	interestVertices <- c(1, 81)
} else if(graph == "grid9Augmented")
{
	graph <- igraph::make_lattice(dim = 2, length = 9)
	graph <- add_vertices(graph, 2)
	for(i in 1:9)
	{
		graph <- add_edges(graph, c(82, i))
		graph <- add_edges(graph, c(83, i+72))
	}
	interestVertices <- c(82, 83)
} else
{
	stop("Unknown graph")
}

if(method == "approxZeroVariance")
{
	counter <- 1
	if(file.exists(outputFile))
	{
		load(outputFile)
		counter <- length(results)+1
	} else results <- list()
	while(counter < nReps + 1)
	{
		results[[counter]] <- approximateZeroVariance(graph = graph, probability = probability, n = sampleSize, seed = counter, interestVertices = interestVertices)
		if(counter %% 100 == 0)
		{
			save(results, file = tmpFile)
			file.rename(from = tmpFile, to = outputFile)
		}
		counter <- counter + 1
	}
} else if(method == "approxZeroVarianceWOR")
{
	counter <- 1
	if(file.exists(outputFile))
	{
		load(outputFile)
		counter <- length(results)+1
	} else results <- list()
	while(counter < nReps + 1)
	{
		results[[counter]] <- approximateZeroVarianceWOR(graph = graph, probability = probability, n = sampleSize, seed = counter, interestVertices = interestVertices)
		if(counter %% 100 == 0)
		{
			save(results, file = tmpFile)
			file.rename(from = tmpFile, to = outputFile)
		}
		counter <- counter + 1
	}
} else
{
	stop("Unrecognized method")
}