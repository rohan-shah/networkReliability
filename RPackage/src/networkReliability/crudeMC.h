#ifndef CRUDEMC_RPACKAGE_HEADER_GUARD
#define CRUDEMC_RPACKAGE_HEADER_GUARD
#include <Rcpp.h>
namespace networkReliability
{
	SEXP crudeMC_igraph(SEXP graph_sexp, SEXP probability_sexp, SEXP n_sexp, SEXP seed_sexp, SEXP interestVertices_sexp);
	SEXP crudeMC_graphAM(SEXP graph_sexp, SEXP probability_sexp, SEXP n_sexp, SEXP seed_sexp, SEXP interestVertices_sexp);
	SEXP crudeMC_graphNEL(SEXP graph_sexp, SEXP probability_sexp, SEXP n_sexp, SEXP seed_sexp, SEXP interestVertices_sexp);
}
#endif
