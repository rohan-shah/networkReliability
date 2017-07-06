#ifndef EXHAUSTIVE_PROBABILITY_HEADER_GUARD
#define EXHAUSTIVE_PROBABILITY_HEADER_GUARD
#include <vector>
#include "includeMPFRNetworkReliability.h"
namespace networkReliability
{
	mpfr_class exhaustiveProbability(std::vector<mpfr_class>& counts, mpfr_class probability, bool countDisconnected);
}
#endif
