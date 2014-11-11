#ifndef TURNIP_HEADER_GUARD
#define TURNIP_HEADER_GUARD
#include <boost/random/mersenne_twister.hpp>
#include <algorithm>
#include "includeMPFR.h"
#include <vector>
#include "Context.h"
namespace networkReliability
{
	struct TurnipEqualRateInput
	{
		TurnipEqualRateInput(boost::mt19937& randomSource, const Context::internalGraph* graph, const std::vector<int>& interestVertices);
		boost::mt19937& randomSource;
		mpfr_class exponentialRate;
		const Context::internalGraph* graph;
		const std::vector<int>& interestVertices;
		std::vector<std::pair<int, int> > edges;
		size_t n;
		mpfr_class estimateFirstMoment, estimateSecondMoment;
		bool warnedStability;
	};
	void turnipEqualRate(TurnipEqualRateInput& input);
}
#endif