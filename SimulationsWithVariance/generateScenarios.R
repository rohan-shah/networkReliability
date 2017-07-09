probabilities <- c(0.99, 0.999, 0.9999)
sampleSize <- c(10, 20)
graphs <- c("dodecahedron", "grid3", "grid4")
scenarios <- expand.grid(method = c("withMerging", "withoutMerging"), sampleSize = sampleSize, stringsAsFactors=FALSE, graph = graphs, nReps = 60000L, probability = probabilities)

scenarios$file <- apply(scenarios, 1, function(x) paste0(x["method"], "-", as.integer(x["sampleSize"]), "-", x["graph"], "-", x["probability"], ".RData", sep=""))
