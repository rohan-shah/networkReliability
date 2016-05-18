#ifndef CONVERT_GRAPH_HEADER_GUARD
#define CONVERT_GRAPH_HEADER_GUARD
#include <Rcpp.h>
#include "context.h"
namespace networkReliability
{
	enum R_GRAPH_TYPE
	{
		IGRAPH, GRAPHNEL, GRAPHAM
	};
	void convertGraph(SEXP graph_sexp, Context::inputGraph& graphRef, R_GRAPH_TYPE graphType);
	void convertGraphIGraph(SEXP graph_sexp, Context::inputGraph& graphRef);
	void convertGraphAM(SEXP graph_sexp, Context::inputGraph& graphRef);
	void convertGraphNEL(SEXP graph_sexp, Context::inputGraph& graphRef);
}
#endif
