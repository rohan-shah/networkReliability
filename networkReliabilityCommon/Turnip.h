#ifndef TURNIP_HEADER_GUARD
#define TURNIP_HEADER_GUARD
#include <boost/random/mersenne_twister.hpp>
#include <algorithm>
#include "includeMPFR.h"
#include <vector>
#include "Context.h"
namespace networkReliability
{
	struct TurnipInput
	{
		TurnipInput(boost::mt19937& randomSource, Context::internalGraph const& graph, const std::vector<int>& interestVertices);
		boost::mt19937& randomSource;
		std::vector<mpfr_class> exponentialRates;
		Context::internalGraph const& graph;
		const std::vector<int>& interestVertices;
		std::vector<std::pair<int, int> > edges;
		size_t n;
		mpfr_class estimateFirstMoment, estimateSecondMoment;
		bool warnedStability;
	};
	void turnip(TurnipInput& input);
}
#endif