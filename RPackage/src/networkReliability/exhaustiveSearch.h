#ifndef EXHAUSTIVE_SEARCH_R_PACKAGE_HEADER_GUARD
#define EXHAUSTIVE_SEARCH_R_PACKAGE_HEADER_GUARD
#include <Rcpp.h>
namespace networkReliability
{
	SEXP exhaustiveSearch_igraph(SEXP graph_sexp, SEXP interestVertices_SEXP, SEXP countDisconnected_sexp);
	SEXP exhaustiveSearch_graphNEL(SEXP graph_sexp, SEXP interestVertices_SEXP, SEXP countDisconnected_sexp);
	SEXP exhaustiveSearch_graphAM(SEXP graph_sexp, SEXP interestVertices_SEXP, SEXP countDisconnected_sexp);
}
#endif
