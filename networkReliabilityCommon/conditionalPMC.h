#ifndef CONDITIONAL_PMC_HEADER_GUARD
#define CONDITIONAL_PMC_HEADER_GUARD
#include <boost/random/mersenne_twister.hpp>
#include <algorithm>
#include "includeMPFR.h"
#include <vector>
#include "Context.h"
namespace networkReliability
{
	struct ConditionalTurnipInput
	{
		ConditionalTurnipInput(boost::mt19937& randomSource, const Context::internalGraph* graph, const std::vector<int>& interestVertices);
		boost::mt19937& randomSource;
		mpfr_class exponentialRate;
		const Context::internalGraph* graph;
		const std::vector<int>& interestVertices;
		std::vector<std::pair<int, int> > edges;
		std::vector<int> edgeCounts;
		//working memory used to select the edges to deactivate
		std::vector<int> workingEdgeCounts;
		//The number of original edges between each connected group, after random edges are fixed as being off.
		std::vector<int> workingEdgeCounts2;
		//exponential rates after edges are removed
		std::vector<mpfr_class> exponentialRates;
		size_t n;
		mpfr_class estimateFirstMoment, estimateSecondMoment;
		bool warnedStability;
		int minimumInoperative;
	};
	void conditionalPMC(ConditionalTurnipInput& input);
}
#endif