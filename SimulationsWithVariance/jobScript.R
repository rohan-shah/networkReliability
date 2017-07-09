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
sampleSize <- scenarios[SCENARIO_INDEX, "sampleSize"]
graph <- scenarios[SCENARIO_INDEX, "graph"]
probability <- scenarios[SCENARIO_INDEX, "probability"]
nReps <- scenarios[SCENARIO_INDEX, "nReps"]
method <- scenarios[SCENARIO_INDEX, "method"]

if(graph == "dodecahedron")
{
	graph <- igraph::read_graph(file = system.file("data/dodecahedron.graphml", package="networkReliability"), format = "graphml")
	interestVertices <- c(1, 20)
} else if(graph == "grid3")
{
	graph <- igraph::make_lattice(dim = 2, length = 3)
	interestVertices <- c(1, 9)
} else if(graph == "grid4")
{
	graph <- igraph::make_lattice(dim = 2, length = 4)
	interestVertices <- c(1, 16)
} else
{
	stop("Unknown graph")
}

counter <- 1
if(file.exists(outputFile))
{
	load(outputFile)
	counter <- length(results)+1
} else results <- list()
if(method == "withMerging")
{
	while(counter < nReps + 1)
	{
		results[[counter]] <- approximateZeroVarianceWORMergeWithVariance(graph = graph, probability = probability, n = sampleSize, seed = counter, interestVertices = interestVertices)
		if(counter %% 100 == 0)
		{
			save(results, file = tmpFile)
			file.rename(from = tmpFile, to = outputFile)
		}
		cat(counter, " / ", nReps, "\n", sep="")
		counter <- counter + 1
	}
} else if(method == "withoutMerging")
{
	while(counter < nReps + 1)
	{
		results[[counter]] <- approximateZeroVarianceWORWithVariance(graph = graph, probability = probability, n = sampleSize, seed = counter, interestVertices = interestVertices)
		if(counter %% 100 == 0)
		{
			save(results, file = tmpFile)
			file.rename(from = tmpFile, to = outputFile)
		}
		cat(counter, " / ", nReps, "\n", sep="")
		counter <- counter + 1
	}
} else
{
	stop("Unrecognized method")
}
save(results, file = tmpFile)
