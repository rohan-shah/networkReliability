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

legendPositions <- vector(mode = "list", length = nrow(scenarios))
for(i in 1:nrow(scenarios))
{
	load(file.path("results", scenarios[i, "file"]))
	endPoints <- seq(step, length(results), by = step)
	allEstimates <- as.numeric(do.call(c, lapply(results, function(x) x@estimate)))
	allVarianceEstimates <- as.numeric(do.call(c, lapply(results, function(x) x@varianceEstimate)))

	empiricalVariances <- sapply(endPoints, function(x) var(allEstimates[1:x]))
	horvitzThomsponVariance <- sapply(endPoints, function(x) mean(allVarianceEstimates[1:x]))
	empiricalMSE <- sapply(endPoints, function(x) mean((allEstimates[1:x] - as.numeric(getExact(scenarios[i, "graph"], scenarios[i, "probability"])))^2))

	data <- rbind(data.frame(x = endPoints, Method = "Empirical Variance", value = empiricalVariances, stringsAsFactors = FALSE), data.frame(x = endPoints, Method = "Horvitz Thompson", value = horvitzThomsponVariance, stringsAsFactors = FALSE))
	#data <- rbind(data, data.frame(x = endPoints, Method = "Empirical MSE", value = empiricalMSE, stringsAsFactors = FALSE))

	pdf(file.path("runningVariances", paste0(scenarios[i, "method"], "-", scenarios[i, "sampleSize"], "-", scenarios[i, "graph"], "-", scenarios[i, "probability"], ".pdf")))
		plot <- ggplot(data = data, mapping = aes(x, value, group = Method, colour = Method)) + geom_line(size = 1) + scale_y_log10(breaks = scales::pretty_breaks(n = 8), labels = scientific_10) + xlab("Number of samples") + theme(axis.text.x = element_text(size = rel(1.5)), axis.text.y = element_text(size = rel(1.25)), axis.title.x = element_text(size = rel(1.5)), axis.title.y = element_blank(), legend.title = element_text(size = rel(1.5)), legend.text = element_text(size = rel(1.5))) + scale_size(guide = "none")
		if(is.null(legendPositions[[i]]))
		{
			plot <- plot + theme(legend.position = c(0.6, 0.2))
		}
		else
		{
		}
		print(plot)
	dev.off()
	cat("Done ", i, " / ", nrow(scenarios), "\n", sep="")
}
