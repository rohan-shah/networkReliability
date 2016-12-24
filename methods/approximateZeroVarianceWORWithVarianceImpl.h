#ifndef APPROXIMATE_ZERO_VARIANCE_WOR_WITH_VARIANCE_HEADER_GUARD
#define APPROXIMATE_ZERO_VARIANCE_WOR_WITH_VARIANCE_HEADER_GUARD
#include "context.h"
namespace networkReliability
{
	struct approximateZeroVarianceWORWithVarianceArgs
	{
	public:
		approximateZeroVarianceWORWithVarianceArgs(const context& contextObj)
			: contextObj(contextObj)
		{}
		const context& contextObj;
		boost::mt19937 randomSource;
		std::size_t n;
		mpfr_class estimate;
		mpfr_class varianceEstimate;
	};
	void approximateZeroVarianceWORWithVariance(approximateZeroVarianceWORWithVarianceArgs& args);
}
#endif
