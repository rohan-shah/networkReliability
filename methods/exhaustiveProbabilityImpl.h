#ifndef EXHAUSTIVE_PROBABILITY_HEADER_GUARD
#define EXHAUSTIVE_PROBABILITY_HEADER_GUARD
#include <vector>
#include "includeMPFRNetworkReliability.h"
namespace networkReliability
{
	double exhaustiveProbability(std::vector<mpfr_class>& counts, mpfr_class probability, bool countDisconnected);
}
#endif
