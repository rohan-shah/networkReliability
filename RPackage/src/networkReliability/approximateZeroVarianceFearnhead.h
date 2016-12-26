#ifndef APPROXIMATE_ZERO_VARIANCE_FEARNHEAD_RPACKAGE_HEADER_GUARD
#define APPROXIMATE_ZERO_VARIANCE_FEARNHEAD_RPACKAGE_HEADER_GUARD
#include <Rcpp.h>
namespace networkReliability
{
	SEXP approximateZeroVarianceFearnhead_igraph(SEXP graph_sexp, SEXP probability_sexp, SEXP n_sexp, SEXP seed_sexp, SEXP interestVertices_sexp);
	SEXP approximateZeroVarianceFearnhead_graphAM(SEXP graph_sexp, SEXP probability_sexp, SEXP n_sexp, SEXP seed_sexp, SEXP interestVertices_sexp);
	SEXP approximateZeroVarianceFearnhead_graphNEL(SEXP graph_sexp, SEXP probability_sexp, SEXP n_sexp, SEXP seed_sexp, SEXP interestVertices_sexp);
}
#endif
