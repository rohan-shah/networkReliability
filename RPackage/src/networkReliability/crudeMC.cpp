#include "crudeMC.h"
#include "crudeMCImpl.h"
#include "convertGraph.h"
#include "createContext.h"
namespace networkReliability
{
	SEXP crudeMC(SEXP graph_sexp, SEXP probability_sexp, SEXP n_sexp, SEXP seed_sexp, SEXP interestVertices_sexp, R_GRAPH_TYPE graphType)
	{
	BEGIN_RCPP
		std::vector<int> interestVertices = Rcpp::as<std::vector<int> >(interestVertices_sexp);
		//Change to base-zero indices
		for(std::vector<int>::iterator i = interestVertices.begin(); i != interestVertices.end(); i++)
		{
			(*i)--;
			if(*i < 0) throw std::runtime_error("Interest vertices must be non-negative");
		}
		double probability = Rcpp::as<double>(probability_sexp);
		mpfr_class probability_mpfr = probability;
		std::size_t n = (std::size_t)Rcpp::as<int>(n_sexp);
		int seed = Rcpp::as<int>(seed_sexp);
		context contextObj = createContext(graph_sexp, interestVertices, probability_mpfr, graphType);

		crudeMCArgs args(contextObj);
		args.randomSource.seed(seed);
		args.n = n;
		return Rcpp::wrap(crudeMC(args));
	END_RCPP
	}
	SEXP crudeMC_igraph(SEXP graph_sexp, SEXP probability_sexp, SEXP n_sexp, SEXP seed_sexp, SEXP interestVertices_sexp)
	{
		return crudeMC(graph_sexp, probability_sexp, n_sexp, seed_sexp, interestVertices_sexp, IGRAPH);
	}
	SEXP crudeMC_graphAM(SEXP graph_sexp, SEXP probability_sexp, SEXP n_sexp, SEXP seed_sexp, SEXP interestVertices_sexp)
	{
		return crudeMC(graph_sexp, probability_sexp, n_sexp, seed_sexp, interestVertices_sexp, GRAPHAM);
	}
	SEXP crudeMC_graphNEL(SEXP graph_sexp, SEXP probability_sexp, SEXP n_sexp, SEXP seed_sexp, SEXP interestVertices_sexp)
	{
		return crudeMC(graph_sexp, probability_sexp, n_sexp, seed_sexp, interestVertices_sexp, GRAPHNEL);
	}
}
