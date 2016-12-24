#ifndef APPROXIMATE_ZERO_VARIANCE_WOR_WITH_VARIANCE_RPACKAGE_HEADER_GUARD
#define APPROXIMATE_ZERO_VARIANCE_WOR_WITH_VARIANCE_RPACKAGE_HEADER_GUARD
#include <Rcpp.h>
namespace networkReliability
{
	SEXP approximateZeroVarianceWORWithVariance_igraph(SEXP graph_sexp, SEXP probability_sexp, SEXP n_sexp, SEXP seed_sexp, SEXP interestVertices_sexp);
	SEXP approximateZeroVarianceWORWithVariance_graphAM(SEXP graph_sexp, SEXP probability_sexp, SEXP n_sexp, SEXP seed_sexp, SEXP interestVertices_sexp);
	SEXP approximateZeroVarianceWORWithVariance_graphNEL(SEXP graph_sexp, SEXP probability_sexp, SEXP n_sexp, SEXP seed_sexp, SEXP interestVertices_sexp);
}
#endif
