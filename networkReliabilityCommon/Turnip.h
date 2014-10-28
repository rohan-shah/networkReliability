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
		TurnipInput(boost::mt19937& randomSource, Context const& context);
		boost::mt19937& randomSource;
		std::vector<mpfr_class> exponentialRates;
		std::vector<std::pair<int, int> > vertices;
		Context const& context;
		size_t n;
		mpfr_class estimateFirstMoment, estimateSecondMoment;
		bool warnedStability;
	};
	void turnip(TurnipInput& input);
}
#endif