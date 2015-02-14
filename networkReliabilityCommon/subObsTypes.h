#ifndef SUB_OBS_TYPE_HEADER_GUARD
#define SUB_OBS_TYPE_HEADER_GUARD
#include "includeMPFR.h"
namespace networkReliability
{
	namespace obs
	{
		struct withImportanceResamplingConstructorType
		{
			int conditioningCount;
			mpfr_class conditioningProb;
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