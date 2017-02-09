#ifndef RESIDUAL_RESAMPLING_PACKAGE_HEADER_GUARD
#define RESIDUAL_RESAMPLING_RPACKAGE_HEADER_GUARD
#include <Rcpp.h>
namespace networkReliability
{
	SEXP residualResampling_igraph(SEXP graph_sexp, SEXP probability_sexp, SEXP n_sexp, SEXP seed_sexp, SEXP interestVertices_sexp);
	SEXP residualResampling_graphAM(SEXP graph_sexp, SEXP probability_sexp, SEXP n_sexp, SEXP seed_sexp, SEXP interestVertices_sexp);
	SEXP residualResampling_graphNEL(SEXP graph_sexp, SEXP probability_sexp, SEXP n_sexp, SEXP seed_sexp, SEXP interestVertices_sexp);
}
#endif
