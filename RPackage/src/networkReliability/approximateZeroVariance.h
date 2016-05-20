#ifndef APPROXIMATE_ZERO_VARIANCE_RPACKAGE_HEADER_GUARD
#define APPROXIMATE_ZERO_VARIANCE_RPACKAGE_HEADER_GUARD
#include <Rcpp.h>
namespace networkReliability
{
	SEXP approximateZeroVariance_igraph(SEXP graph_sexp, SEXP probability_sexp, SEXP n_sexp, SEXP seed_sexp, SEXP interestVertices_sexp);
	SEXP approximateZeroVariance_graphAM(SEXP graph_sexp, SEXP probability_sexp, SEXP n_sexp, SEXP seed_sexp, SEXP interestVertices_sexp);
	SEXP approximateZeroVariance_graphNEL(SEXP graph_sexp, SEXP probability_sexp, SEXP n_sexp, SEXP seed_sexp, SEXP interestVertices_sexp);
}
#endif
