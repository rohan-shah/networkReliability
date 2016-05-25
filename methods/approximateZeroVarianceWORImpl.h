#ifndef APPROXIMATE_ZERO_VARIANCE_WOR_HEADER_GUARD
#define APPROXIMATE_ZERO_VARIANCE_WOR_HEADER_GUARD
#include "context.h"
namespace networkReliability
{
	struct approximateZeroVarianceWORArgs
	{
	public:
		approximateZeroVarianceWORArgs(const context& contextObj)
			: contextObj(contextObj)
		{}
		const context& contextObj;
		boost::mt19937 randomSource;
		std::size_t n;
		mpfr_class estimate;
		bool optimiseMinCut;
	};
	void approximateZeroVarianceWOR(approximateZeroVarianceWORArgs& args);
}
#endif
