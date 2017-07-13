#ifndef APPROXIMATE_ZERO_VARIANCE_WOR_MERGE_WITH_VARIANCE_HEADER_GUARD
#define APPROXIMATE_ZERO_VARIANCE_WOR_MERGE_WITH_VARIANCE_HEADER_GUARD
#include "context.h"
#include "approximateZeroVarianceWORMergeImpl.h"
#include "approximateZeroVarianceWORImpl.h"
namespace networkReliability
{
	struct approximateZeroVarianceWORMergeWithVarianceArgs
	{
	public:
		approximateZeroVarianceWORMergeWithVarianceArgs(const context& contextObj)
			: contextObj(contextObj), graphFile("")
		{}
		const context& contextObj;
		boost::mt19937 randomSource;
		std::size_t n;
		mpfr_class estimate, varianceEstimate;
		std::string graphFile;
	};
	void approximateZeroVarianceWORMergeWithVariance(approximateZeroVarianceWORMergeWithVarianceArgs& args);
}
#endif
