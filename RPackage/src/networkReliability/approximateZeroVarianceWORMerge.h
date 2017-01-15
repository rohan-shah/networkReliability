#ifndef APPROXIMATE_ZERO_VARIANCE_WOR_MERGE_RPACKAGE_HEADER_GUARD
#define APPROXIMATE_ZERO_VARIANCE_WOR_MERGE_RPACKAGE_HEADER_GUARD
#include <Rcpp.h>
namespace networkReliability
{
	SEXP approximateZeroVarianceWORMerge_igraph(SEXP graph_sexp, SEXP probability_sexp, SEXP n_sexp, SEXP seed_sexp, SEXP interestVertices_sexp);
	SEXP approximateZeroVarianceWORMerge_graphAM(SEXP graph_sexp, SEXP probability_sexp, SEXP n_sexp, SEXP seed_sexp, SEXP interestVertices_sexp);
	SEXP approximateZeroVarianceWORMerge_graphNEL(SEXP graph_sexp, SEXP probability_sexp, SEXP n_sexp, SEXP seed_sexp, SEXP interestVertices_sexp);
}
#endif
