#ifndef RESIDUAL_RESAMPLING_HEADER_GUARD
#define RESIDUAL_RESAMPLING_HEADER_GUARD
#include "context.h"
namespace networkReliability
{
	struct residualResamplingArgs
	{
	public:
		residualResamplingArgs(const context& contextObj)
			: contextObj(contextObj)
		{}
		const context& contextObj;
		boost::mt19937 randomSource;
		std::size_t n;
		mpfr_class estimate;
	};
	void residualResampling(residualResamplingArgs& args);
}
#endif
