load("./summarised.RData")
source("./generateScenarios.R")
library(ggplot2)
library(dplyr)
size <- 1
data <- cbind(scenarios, wnrv = wnrv, relativeError = relativeErrors, empiricalVariance = empiricalVariances, workNormalizedVariance = workNormalizedVariance)
scale_y <- scale_y_log10(label = trans_format("log10", math_format(10^.x)))

#Plot workNormalizedVariance, for the dodecahedron with p = 0.99
plot <- data %>% filter(graph == "dodecahedron" & probability == 0.99 & method != "approxZeroVarianceWORWithVariance") %>% ggplot(mapping = aes(sampleSize, workNormalizedVariance, colour = method)) + scale_y + geom_line(size = size) + scale_x_log10(breaks = c(10, 100, 1000, 10000), labels = c(10, 100, 1000, 10000)) + theme(legend.position = c(0.2, 0.2)) + xlab("Sample size") + ylab("WNV") + scale_colour_discrete("Method", breaks = c("approxZeroVariance", "approxZeroVarianceWOR", "approxZeroVarianceWORMerge", "fearnhead"), labels = c("IS", "WOR", "WOR-Merge", "Fearnhead"))
pdf("./dodecahedronWNV0_99.pdf")
print(plot)
dev.off()
#Plot RE, for the dodecahedron with p = 0.99
plot <- data %>% filter(graph == "dodecahedron" & probability == 0.99 & method != "approxZeroVarianceWORWithVariance") %>% ggplot(mapping = aes(sampleSize, relativeError, colour = method)) + scale_y + geom_line(size = size) + scale_x_log10(breaks = c(10, 100, 1000, 10000), labels = c(10, 100, 1000, 10000)) + theme(legend.position = c(0.2, 0.2)) + xlab("Sample size") + ylab("RE") + scale_colour_discrete("Method", breaks = c("approxZeroVariance", "approxZeroVarianceWOR", "approxZeroVarianceWORMerge", "fearnhead"), labels = c("IS", "WOR", "WOR-Merge", "Fearnhead"))
pdf("./dodecahedronRE0_99.pdf")
print(plot)
dev.off()


#Plot workNormalizedVariance, for the augmented grid graph with p = 0.99
plot <- data %>% filter(graph == "grid9Augmented" & probability == 0.99 & method != "approxZeroVarianceWORWithVariance") %>% ggplot(mapping = aes(sampleSize, workNormalizedVariance, colour = method)) + scale_y + geom_line(size = size) + scale_x_log10(breaks = c(10, 100, 1000, 10000), labels = c(10, 100, 1000, 10000)) + theme(legend.position = c(0.2, 0.2)) + xlab("Sample size") + ylab("WNV") + scale_colour_discrete("Method", breaks = c("approxZeroVariance", "approxZeroVarianceWOR", "approxZeroVarianceWORMerge", "fearnhead"), labels = c("IS", "WOR", "WOR-Merge", "Fearnhead"))
pdf("./grid9AugmentedWNV0_99.pdf")
print(plot)
dev.off()
#Plot RE, for the augmented grid graph with p = 0.99
plot <- data %>% filter(graph == "grid9Augmented" & probability == 0.99 & method != "approxZeroVarianceWORWithVariance") %>% ggplot(mapping = aes(sampleSize, relativeError, colour = method)) + scale_y + geom_line(size = size) + scale_x_log10(breaks = c(10, 100, 1000, 10000), labels = c(10, 100, 1000, 10000)) + theme(legend.position = c(0.2, 0.2)) + xlab("Sample size") + ylab("RE") + scale_colour_discrete("Method", breaks = c("approxZeroVariance", "approxZeroVarianceWOR", "approxZeroVarianceWORMerge", "fearnhead"), labels = c("IS", "WOR", "WOR-Merge", "Fearnhead"))
pdf("./grid9AugmentedRE0_99.pdf")
print(plot)
dev.off()

#Plot workNormalizedVariance, for the grid graph with p = 0.99
plot <- data %>% filter(graph == "grid9" & probability == 0.99 & method != "approxZeroVarianceWORWithVariance") %>% ggplot(mapping = aes(sampleSize, workNormalizedVariance, colour = method)) + scale_y + geom_line(size = size) + scale_x_log10(breaks = c(10, 100, 1000, 10000), labels = c(10, 100, 1000, 10000)) + theme(legend.position = c(0.2, 0.2)) + xlab("Sample size") + ylab("WNV") + scale_colour_discrete("Method", breaks = c("approxZeroVariance", "approxZeroVarianceWOR", "approxZeroVarianceWORMerge", "fearnhead"), labels = c("IS", "WOR", "WOR-Merge", "Fearnhead"))
pdf("./grid9WNV0_99.pdf")
print(plot)
dev.off()
#Plot RE, for the grid graph with p = 0.99
plot <- data %>% filter(graph == "grid9" & probability == 0.99 & method != "approxZeroVarianceWORWithVariance") %>% ggplot(mapping = aes(sampleSize, relativeError, colour = method)) + scale_y + geom_line(size = size) + scale_x_log10(breaks = c(10, 100, 1000, 10000), labels = c(10, 100, 1000, 10000)) + theme(legend.position = c(0.2, 0.2)) + xlab("Sample size") + ylab("RE") + scale_colour_discrete("Method", breaks = c("approxZeroVariance", "approxZeroVarianceWOR", "approxZeroVarianceWORMerge", "fearnhead"), labels = c("IS", "WOR", "WOR-Merge", "Fearnhead"))
pdf("./grid9RE0_99.pdf")
print(plot)
dev.off()

#Plot workNormalizedVariance, for the dodecahedron series graph with p = 0.9999
plot <- data %>% filter(graph == "dodecSeries" & probability == 0.9999 & method != "approxZeroVarianceWORWithVariance") %>% ggplot(mapping = aes(sampleSize, workNormalizedVariance, colour = method)) + scale_y + geom_line(size = size) + scale_x_log10(breaks = c(10, 100, 1000, 10000), labels = c(10, 100, 1000, 10000)) + theme(legend.position = c(0.2, 0.2)) + xlab("Sample size") + ylab("WNV") + scale_colour_discrete("Method", breaks = c("approxZeroVariance", "approxZeroVarianceWOR", "approxZeroVarianceWORMerge", "fearnhead"), labels = c("IS", "WOR", "WOR-Merge", "Fearnhead"))
pdf("./dodecSeriesWNV0_9999.pdf")
print(plot)
dev.off()
#Plot RE, for the dodecahedron series graph with p = 0.9999
plot <- data %>% filter(graph == "dodecSeries" & probability == 0.9999 & method != "approxZeroVarianceWORWithVariance") %>% ggplot(mapping = aes(sampleSize, relativeError, colour = method)) + scale_y + geom_line(size = size) + scale_x_log10(breaks = c(10, 100, 1000, 10000), labels = c(10, 100, 1000, 10000)) + theme(legend.position = c(0.2, 0.2)) + xlab("Sample size") + ylab("RE") + scale_colour_discrete("Method", breaks = c("approxZeroVariance", "approxZeroVarianceWOR", "approxZeroVarianceWORMerge", "fearnhead"), labels = c("IS", "WOR", "WOR-Merge", "Fearnhead"))
pdf("./dodecSeriesRE0_9999.pdf")
print(plot)
dev.off()
#Plot variance, for the dodecahedron series graph with p = 0.9999
plot <- data %>% filter(graph == "dodecSeries" & probability == 0.9999 & method != "approxZeroVarianceWORWithVariance") %>% ggplot(mapping = aes(sampleSize, empiricalVariance, colour = method)) + scale_y + geom_line(size = size) + scale_x_log10(breaks = c(10, 100, 1000, 10000), labels = c(10, 100, 1000, 10000)) + theme(legend.position = c(0.2, 0.2)) + xlab("Sample size") + ylab("Var") + scale_colour_discrete("Method", breaks = c("approxZeroVariance", "approxZeroVarianceWOR", "approxZeroVarianceWORMerge", "fearnhead"), labels = c("IS", "WOR", "WOR-Merge", "Fearnhead"))
pdf("./dodecSeriesVar0_9999.pdf")
print(plot)
dev.off()

#Plot workNormalizedVariance, for the dodecahedron parallel graph with p = 0.9999
plot <- data %>% filter(graph == "dodecParallel" & probability == 0.9999 & method != "approxZeroVarianceWORWithVariance") %>% ggplot(mapping = aes(sampleSize, workNormalizedVariance, colour = method)) + scale_y + geom_line(size = size) + scale_x_log10(breaks = c(10, 100, 1000, 10000), labels = c(10, 100, 1000, 10000)) + theme(legend.position = c(0.2, 0.2)) + xlab("Sample size") + ylab("WNV") + scale_colour_discrete("Method", breaks = c("approxZeroVariance", "approxZeroVarianceWOR", "approxZeroVarianceWORMerge", "fearnhead"), labels = c("IS", "WOR", "WOR-Merge", "Fearnhead"))
pdf("./dodecParallelWNV0_9999.pdf")
print(plot)
dev.off()
#Plot RE, for the dodecahedron parallel graph with p = 0.9999
plot <- data %>% filter(graph == "dodecParallel" & probability == 0.9999 & method != "approxZeroVarianceWORWithVariance") %>% ggplot(mapping = aes(sampleSize, relativeError, colour = method)) + scale_y + geom_line(size = size) + scale_x_log10(breaks = c(10, 100, 1000, 10000), labels = c(10, 100, 1000, 10000)) + theme(legend.position = c(0.2, 0.2)) + xlab("Sample size") + ylab("RE") + scale_colour_discrete("Method", breaks = c("approxZeroVariance", "approxZeroVarianceWOR", "approxZeroVarianceWORMerge", "fearnhead"), labels = c("IS", "WOR", "WOR-Merge", "Fearnhead"))
pdf("./dodecParallelRE0_9999.pdf")
print(plot)
dev.off()