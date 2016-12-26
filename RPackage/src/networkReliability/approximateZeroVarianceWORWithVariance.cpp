#include "approximateZeroVarianceWORWithVariance.h"
#include "approximateZeroVarianceWORWithVarianceImpl.h"
#include "convertGraph.h"
#include "createContext.h"
namespace networkReliability
{
	SEXP approximateZeroVarianceWORWithVariance(SEXP graph_sexp, SEXP probability_sexp, SEXP n_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, R_GRAPH_TYPE graphType)
	{
	BEGIN_RCPP
		std::vector<int> interestVertices = Rcpp::as<std::vector<int> >(interestVertices_sexp);
		//Change to base-zero indices
		for(std::vector<int>::iterator i = interestVertices.begin(); i != interestVertices.end(); i++)
		{
			(*i)--;
		}
		double probability = Rcpp::as<double>(probability_sexp);
		mpfr_class probability_mpfr = probability;
		std::size_t n = (std::size_t)Rcpp::as<int>(n_sexp);
		int seed = Rcpp::as<int>(seed_sexp);
		context contextObj = createContext(graph_sexp, interestVertices, probability_mpfr, graphType);

		approximateZeroVarianceWORWithVarianceArgs args(contextObj);
		args.randomSource.seed(seed);
		args.n = n;
		approximateZeroVarianceWORWithVariance(args);
		return Rcpp::List::create(Rcpp::Named("estimate") = Rcpp::wrap(args.estimate.str()), Rcpp::Named("varianceEstimate") = Rcpp::wrap(args.varianceEstimate.str()));
	END_RCPP
	}
	SEXP approximateZeroVarianceWORWithVariance_igraph(SEXP graph_sexp, SEXP probability_sexp, SEXP n_sexp, SEXP seed_sexp, SEXP interestVertices_sexp)
	{
		return approximateZeroVarianceWORWithVariance(graph_sexp, probability_sexp, n_sexp, seed_sexp, interestVertices_sexp, IGRAPH);
	}
	SEXP approximateZeroVarianceWORWithVariance_graphAM(SEXP graph_sexp, SEXP probability_sexp, SEXP n_sexp, SEXP seed_sexp, SEXP interestVertices_sexp)
	{
		return approximateZeroVarianceWORWithVariance(graph_sexp, probability_sexp, n_sexp, seed_sexp, interestVertices_sexp, GRAPHAM);
	}
	SEXP approximateZeroVarianceWORWithVariance_graphNEL(SEXP graph_sexp, SEXP probability_sexp, SEXP n_sexp, SEXP seed_sexp, SEXP interestVertices_sexp)
	{
		return approximateZeroVarianceWORWithVariance(graph_sexp, probability_sexp, n_sexp, seed_sexp, interestVertices_sexp, GRAPHNEL);
	}
}
