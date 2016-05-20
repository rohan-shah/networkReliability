#ifndef CREATE_CONTEXT_HEADER_GUARD
#define CREATE_CONTEXT_HEADER_GUARD
#include "Rcpp.h"
#include "context.h"
#include "convertGraph.h"
namespace networkReliability
{
	context createContext(SEXP graph_sexp, std::vector<int>& interestVertices, mpfr_class probability, R_GRAPH_TYPE type);
}
#endif
