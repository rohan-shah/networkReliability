#ifndef TURNIP_HEADER_GUARD
#define TURNIP_HEADER_GUARD
#include <boost/random/mersenne_twister.hpp>
#include <algorithm>
#include "includeMPFRNetworkReliability.h"
#include <vector>
#include "context.h"
namespace networkReliability
{
	struct turnipEqualRateInput
	{
		turnipEqualRateInput(boost::mt19937& randomSource, const context::internalGraph* graph, const std::vector<int>& interestVertices);
		boost::mt19937& randomSource;
		mpfr_class exponentialRate;
		const context::internalGraph* graph;
		const std::vector<int>& interestVertices;
		std::vector<std::pair<int, int> > edges;
		size_t n;
		mpfr_class estimateFirstMoment, estimateSecondMoment;
		bool warnedStability;
	};
	void turnipEqualRate(turnipEqualRateInput& input);
}
#endif
