library(Rmpfr)
hasPackage <- require(networkReliability)
library(R.utils)
if(!hasPackage) q(save="no")
if(!file.exists("data/dodecahedronConnectedEnumeration.RData"))
{
	dodecahedron <- igraph::read_graph(file = system.file("data/dodecahedron.graphml", package="networkReliability"), format = "graphml")
	dodecahedronConnectedEnumeration <- exhaustiveSearch(graph = dodecahedron, interestVertices = c(1, 20), countDisconnected = FALSE)
	save(dodecahedronConnectedEnumeration, file = "data/dodecahedronConnectedEnumeration.RData")
}
if(!file.exists("data/dodecahedronDisconnectedEnumeration.RData"))
{
	dodecahedron <- igraph::read_graph(file = system.file("data/dodecahedron.graphml", package="networkReliability"), format = "graphml")
	dodecahedronDisconnectedEnumeration <- exhaustiveSearch(graph = dodecahedron, interestVertices = c(1, 20), countDisconnected = TRUE)
	save(dodecahedronDisconnectedEnumeration, file = "data/dodecahedronDisconnectedEnumeration.RData")
}
