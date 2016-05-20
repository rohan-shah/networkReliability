#ifndef CRUDEMC_IMPL_HEADER_GUARD
#define CRUDEMC_IMPL_HEADER_GUARD
#include "context.h"
namespace networkReliability
{
	struct crudeMCArgs
	{
	public:
		crudeMCArgs(const context& contextObj)
			:contextObj(contextObj)
		{}
		const context& contextObj;
		boost::mt19937 randomSource;
		std::size_t n;
	};
	double crudeMC(crudeMCArgs& args);
}
#endif
