#include "computeConditionalProb.h"
#include "mpfr.h"
namespace networkReliability
{
	mpfr_class computeConditionalProb(const std::vector<mpfr_class>& rates)
	{
		mpfr_class ret = 0;
		for(std::vector<mpfr_class>::const_iterator i = rates.begin(); i != rates.end(); i++)
		{
			mpfr_class part = mpfr::exp(-(*i), MPFR_RNDN);
			for(std::vector<mpfr_class>::const_iterator j = rates.begin(); j != rates.end(); j++)
			{
				if(j != i)
				{
					part *= *j / (*j - *i);
				}
			}
			ret += part;
		}
		return ret;
	}
	mpfr_class computeConditionalProb(const std::vector<mpfr_class>& rates, std::vector<mpfr_class>& scratchMemory)
	{
		mpfr_class productRates = 1;
		for (std::vector<mpfr_class>::const_iterator i = rates.begin(); i != rates.end(); i++)
		{
			productRates *= *i;
		}
		scratchMemory.resize(rates.size());
		for (std::vector<mpfr_class>::const_iterator i = rates.begin(); i != rates.end(); i++)
		{
			scratchMemory[std::distance(rates.begin(), i)] = mpfr::exp(-(*i), MPFR_RNDN) * productRates / *i;
		}
		for (std::vector<mpfr_class>::const_iterator i = rates.begin(); i != rates.end(); i++)
		{
			for (std::vector<mpfr_class>::const_iterator j = i + 1; j != rates.end(); j++)
			{
				mpfr_class inverse = 1 / (*j - *i);
				scratchMemory[std::distance(rates.begin(), i)] *= inverse;
				scratchMemory[std::distance(rates.begin(), j)] *= -inverse;
			}
		}
		mpfr_class ret = 0;
		for (std::vector<mpfr_class>::const_iterator i = scratchMemory.begin(); i != scratchMemory.end(); i++)
		{
			ret += *i;
		}
		return ret;
	}
}