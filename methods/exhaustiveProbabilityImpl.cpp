#include "exhaustiveProbabilityImpl.h"
#include <boost/filesystem.hpp>
namespace networkReliability
{
	mpfr_class exhaustiveProbability(std::vector<mpfr_class>& counts, mpfr_class probability, bool countDisconnected)
	{
		if(probability < 0 || probability > 1)
		{
			throw std::runtime_error("Input probability must be a value in (0, 1)");
		}
		mpfr_class compProbability = 1 - probability;
		const std::size_t nEdges = counts.size() - 1;

		mpfr_class result = 0;
		for(std::size_t i = 0; i < nEdges+1; i++)
		{
			mpfr_class probabilityPower = boost::multiprecision::pow(probability, i);
			mpfr_class compProbabilityPower = boost::multiprecision::pow(compProbability, (unsigned long)(nEdges - i));

			result += counts[i] * compProbabilityPower * probabilityPower;
		}
		mpfr_class unreliability;
		if(countDisconnected)
		{
			unreliability = result;
		}
		else
		{
			unreliability = 1 - result;
		}
		return unreliability;
	}
}
