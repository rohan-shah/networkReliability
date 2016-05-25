#include "setDefaultPrec.h"
#include "includeMPFRNetworkReliability.h"
namespace networkReliability
{
	SEXP setDefaultPrec(SEXP prec_sexp)
	{
	BEGIN_RCPP
		networkReliability::mpfr_class::default_precision(Rcpp::as<int>(prec_sexp));
	VOID_END_RCPP
		return R_NilValue;
	}
}
