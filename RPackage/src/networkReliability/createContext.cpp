#include "createContext.h"
namespace networkReliability
{
	context createContext(SEXP graph_sexp, std::vector<int>& interestVertices, mpfr_class probability, R_GRAPH_TYPE type)
	{
		boost::shared_ptr<context::inputGraph> graph(new context::inputGraph());
		convertGraph(graph_sexp, *graph, type);

		std::size_t nEdges = boost::num_edges(*graph);
		boost::shared_ptr<std::vector<unsigned int> > edgeOrdering(new std::vector<unsigned int>(nEdges));
		for(int i = 0; i < (int)nEdges; i++) (*edgeOrdering)[i] = (unsigned int)i;

		boost::shared_ptr<std::vector<int> > interestVertices_copy(new std::vector<int>());
		interestVertices_copy->insert(interestVertices_copy->begin(), interestVertices.begin(), interestVertices.end());

		boost::shared_ptr<std::vector<context::vertexPosition> > vertexPositions;
		return context(graph, edgeOrdering, interestVertices_copy, vertexPositions, probability);
	}
}
