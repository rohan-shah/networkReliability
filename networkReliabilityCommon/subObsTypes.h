#ifndef SUB_OBS_TYPE_HEADER_GUARD
#define SUB_OBS_TYPE_HEADER_GUARD
#include "includeMPFRNetworkReliability.h"
namespace networkReliability
{
	namespace obs
	{
		struct withImportanceResamplingConstructorType
		{
			int conditioningCount;
			mpfr_class conditioningProb;
			mpfr_class resamplingProb;
		};
		struct withResamplingConstructorType
		{
			int conditioningCount;
			mpfr_class conditioningProb;
		};
		struct basicConstructorType
		{
		};
	}
	namespace subObs
	{
		struct withImportanceResamplingConstructorType
		{
			int conditioningCount;
			mpfr_class conditioningProb;
			mpfr_class resamplingProb;
			double nextRadius;
			std::vector<int> boundaryEdges;
			std::vector<int> interiorEdges;
		};
		struct withResamplingConstructorType
		{
			int conditioningCount;
			mpfr_class conditioningProb;
		};
		struct basicConstructorType
		{
		};
	}
}
#endif
