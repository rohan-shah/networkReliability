methods <- c("approxZeroVariance", "approxZeroVarianceWOR")
probabilities <- c(0.9, 0.99, 0.999, 0.9999)
sampleSize <- c(10, 100, 10000)
graphs <- c("dodecahedron", "grid9", "grid9Augmented", "dodecSeries", "dodecParallel")
scenarios <- expand.grid(method = methods, sampleSize = sampleSize, stringsAsFactors=FALSE, graph = graphs, nReps = 5000L, probability = probabilities)

scenarios$file <- apply(scenarios, 1, function(x) paste0(x["method"], "-", as.integer(x["sampleSize"]), "-", x["graph"], "-", x["probability"], ".RData", sep=""))
