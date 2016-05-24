library(Rmpfr)
hasPackage <- require(networkReliability)
library(R.utils)
if(!hasPackage) q(save="no")
if(!file.exists("data/dodecahedronEnumeration.RData"))
{
	dodecahedron <- igraph::read_graph(file = system.file("data/dodecahedron.graphml", package="networkReliability"), format = "graphml")
	dodecahedronEnumeration <- exhaustiveSearch(graph = dodecahedron, interestVertices = c(1, 20), countDisconnected = FALSE)
	save(dodecahedronEnumeration, file = "data/dodecahedronEnumeration.RData")
}
