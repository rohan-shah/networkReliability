#include "exhaustiveSearch.h"
#include "exhaustiveSearchImpl.h"
#include "convertGraph.h"
namespace networkReliability
{
	SEXP exhaustiveSearch(SEXP graph_sexp, SEXP interestVertices_sexp, SEXP countDisconnected_sexp, R_GRAPH_TYPE graphType)
	{
	BEGIN_RCPP
		std::vector<int> interestVertices = Rcpp::as<std::vector<int> >(interestVertices_sexp);
		context::inputGraph graph;
		convertGraph(graph_sexp, graph, graphType);
		std::size_t nVertices = boost::num_vertices(graph);
		//Change to base-zero indices
		for(std::vector<int>::iterator i = interestVertices.begin(); i != interestVertices.end(); i++)
		{
			(*i)--;
			if(*i < 0) throw std::runtime_error("Interest vertices must be non-negative");
			if(*i >= (int)nVertices) throw std::runtime_error("Input vertex was too large");
		}

		exhaustiveSearchArgs args(graph);
		args.interestVertices = interestVertices;
		args.countDisconnected = Rcpp::as<bool>(countDisconnected_sexp);
		exhaustiveSearch(args);

		Rcpp::CharacterVector result(args.result.size());
		for(std::vector<exhaustiveSearchArgs::counterType>::iterator i = args.result.begin(); i != args.result.end(); i++)
		{
			std::stringstream ss;
			ss << *i;
			result[std::distance(args.result.begin(), i)] = ss.str();
		}
		return result;
	END_RCPP
	}
	SEXP exhaustiveSearch_igraph(SEXP graph_sexp, SEXP interestVertices_sexp, SEXP countDisconnected_sexp)
	{
		return exhaustiveSearch(graph_sexp, interestVertices_sexp, countDisconnected_sexp, IGRAPH);
	}
	SEXP exhaustiveSearch_graphNEL(SEXP graph_sexp, SEXP interestVertices_sexp, SEXP countDisconnected_sexp)
	{
		return exhaustiveSearch(graph_sexp, interestVertices_sexp, countDisconnected_sexp, GRAPHNEL);
	}
	SEXP exhaustiveSearch_graphAM(SEXP graph_sexp, SEXP interestVertices_sexp, SEXP countDisconnected_sexp)
	{
		return exhaustiveSearch(graph_sexp, interestVertices_sexp, countDisconnected_sexp, GRAPHAM);
	}
}
