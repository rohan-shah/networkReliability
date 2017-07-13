#ifndef APPROXIMATE_ZERO_VARIANCE_HEADER_GUARD
#define APPROXIMATE_ZERO_VARIANCE_HEADER_GUARD
#include "context.h"
#include <boost/random/mersenne_twister.hpp>
namespace networkReliability
{
	struct approximateZeroVarianceArgs
	{
	public:
		approximateZeroVarianceArgs(const context& contextObj)
			: contextObj(contextObj)
		{}
		const context& contextObj;
		boost::mt19937 randomSource;
		std::size_t n;
		mpfr_class estimateFirstMoment, estimateSecondMoment;
		bool optimiseMinCut;
	};
	void approximateZeroVariance(approximateZeroVarianceArgs& args);
}
#endif
