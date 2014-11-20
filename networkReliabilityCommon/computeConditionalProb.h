#ifndef COMPUTE_CONDITIONAL_PROB_HEADER_GUARD
#define COMPUTE_CONDITIONAL_PROB_HEADER_GUARD
#include <algorithm>
#include "includeMPFR.h"
#include <vector>
namespace networkReliability
{
	mpfr_class computeConditionalProb(const std::vector<mpfr_class>& rates);
	mpfr_class computeConditionalProb(const std::vector<mpfr_class>& rates, std::vector<mpfr_class>& scratchMemory);
	mpfr_class computeConditionalProbLarger(const std::vector<mpfr_class>& firstRates, const std::vector<mpfr_class>& secondRates, std::vector<mpfr_class>& scratchMemory);
}
#endif