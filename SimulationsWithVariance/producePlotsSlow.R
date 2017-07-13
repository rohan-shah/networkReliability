source("./generateScenarios.R")
library(networkReliability)
library(Rmpfr)
library(ggplot2)
step <- 100
if(!dir.exists("runningVariances")) dir.create("runningVariances", showWarnings = FALSE)

data(dodecahedronConnectedEnumeration, package="networkReliability")
if(!exists("grid3Enumeration"))
{
	graph <- igraph::make_lattice(dim = 2, length = 3)
	grid3Enumeration <- networkReliability::exhaustiveSearch(graph = graph, interestVertices = c(1, 9), countDisconnected = FALSE)
}
if(!exists("grid4Enumeration"))
{
	graph <- igraph::make_lattice(dim = 2, length = 4)
	grid4Enumeration <- networkReliability::exhaustiveSearch(graph = graph, interestVertices = c(1, 16), countDisconnected = FALSE)
}
getExact <- function(graph, probability)
{
	if(graph == "dodecahedron")
	{
		return(networkReliability::exhaustiveProbability(dodecahedronConnectedEnumeration, probability))
	} 
	else if(graph == "grid3")
	{
		return(networkReliability::exhaustiveProbability(grid3Enumeration, probability))
	} 
	else if(graph == "grid4")
	{
		return(networkReliability::exhaustiveProbability(grid4Enumeration, probability))
	}
	else
	{
		stop("Unknown graph")
	}
}

scientific_10 <- function(x) 
{
	parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}
for(i in 1:nrow(scenarios))
{
	load(file.path("results", scenarios[i, "file"]))
	endPoints <- seq(step, length(results), by = step)
	allEstimates <- do.call(c, lapply(results, function(x) x@estimate))
	allEstimatesSquared <- allEstimates^2
	allVarianceEstimates <- do.call(c, lapply(results, function(x) x@varianceEstimate))

	empiricalVariances <- sapply(endPoints, function(x) as.numeric(mean(allEstimatesSquared[1:x]) - mean(allEstimates[1:x])^2))
	horvitzThomsponVariance <- sapply(endPoints, function(x) as.numeric(mean(allVarianceEstimates[1:x])))
	empiricalMSE <- sapply(endPoints, function(x)
		{
			firstMoment <- mean(allEstimates[1:x])
			exact <- getExact(scenarios[i, "graph"], scenarios[i, "probability"])
			as.numeric(mean(allEstimatesSquared[1:x]) - 2 * firstMoment * exact + exact^2)
		})

	data <- rbind(data.frame(x = endPoints, method = "Empirical Variance", value = empiricalVariances, stringsAsFactors = FALSE), data.frame(x = endPoints, method = "Horvitz Thompson", value = horvitzThomsponVariance, stringsAsFactors = FALSE), data.frame(x = endPoints, method = "Empirical MSE", value = empiricalMSE, stringsAsFactors = FALSE))

	pdf(file.path("runningVariances", paste0(scenarios[i, "method"], "-", scenarios[i, "sampleSize"], "-", scenarios[i, "graph"], "-", scenarios[i, "probability"], ".pdf")))
		print(ggplot(data = data, mapping = aes(x, value, group = method, colour = method)) + geom_line() + scale_y_log10(breaks = scales::pretty_breaks(n = 8), labels = scientific_10) + ylab("Log10 of estimated variance") + xlab("Number of samples"))
	dev.off()
	cat("Done ", i, " / ", nrow(scenarios), "\n", sep="")
}
