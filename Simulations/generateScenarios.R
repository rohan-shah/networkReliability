methods <- c("approxZeroVariance", "approxZeroVarianceWOR", "fearnhead", "approxZeroVarianceWORMerge")
probabilities <- c(0.99, 0.9999)
sampleSize <- c(10, 20, 100, 1000, 10000)
graphs <- c("dodecahedron", "grid9", "grid9Augmented", "dodecSeries", "dodecParallel")
scenarios <- expand.grid(method = methods, sampleSize = sampleSize, stringsAsFactors=FALSE, graph = graphs, nReps = 1000L, probability = probabilities)

withVarianceScenarios <- expand.grid(method = "approxZeroVarianceWORWithVariance", sampleSize = c(10, 20), stringsAsFactors=FALSE, graph = graphs, nReps = 1000L, probability = probabilities)
scenarios <- rbind(scenarios, withVarianceScenarios)

scenarios$file <- apply(scenarios, 1, function(x) paste0(x["method"], "-", as.integer(x["sampleSize"]), "-", x["graph"], "-", x["probability"], ".RData", sep=""))
