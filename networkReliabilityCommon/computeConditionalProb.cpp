#include "computeConditionalProb.h"
#include "mpfr.h"
namespace networkReliability
{
	mpfr_class computeConditionalProb(const std::vector<mpfr_class>& rates)
	{
		mpfr_class ret = 0;
		for(std::vector<mpfr_class>::const_iterator i = rates.begin(); i != rates.end(); i++)
		{
			mpfr_class part = boost::multiprecision::exp(-(*i));
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
			scratchMemory[std::distance(rates.begin(), i)] = boost::multiprecision::exp(-(*i)) * productRates / *i;
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
	mpfr_class computeConditionalProbLarger(const std::vector<mpfr_class>& firstRates, const std::vector<mpfr_class>& secondRates, std::vector<mpfr_class>& scratchMemory)
	{
		mpfr_class productRates = 1;
		for (std::vector<mpfr_class>::const_iterator i = firstRates.begin(); i != firstRates.end(); i++)
		{
			productRates *= *i;
		}
		for (std::vector<mpfr_class>::const_iterator i = secondRates.begin(); i != secondRates.end(); i++)
		{
			productRates *= *i;
		}
		scratchMemory.resize(firstRates.size() + secondRates.size());

		for (std::vector<mpfr_class>::const_iterator i = secondRates.begin(); i != secondRates.end(); i++)
		{
			scratchMemory[std::distance(secondRates.begin(), i) + firstRates.size()] = boost::multiprecision::exp(-(*i)) * productRates / *i;
		}
		for (std::vector<mpfr_class>::const_iterator i = secondRates.begin(); i != secondRates.end(); i++)
		{
			//fix up numerator of product of rates in previous computation
			for (std::vector<mpfr_class>::iterator scratchIterator = scratchMemory.begin(), endIterator = scratchMemory.begin() + firstRates.size(); scratchIterator != endIterator; scratchIterator++)
			{
				*scratchIterator *= *i;
			}
		}
		//fix up denominator of product of rates in previous computation
		for (std::vector<mpfr_class>::const_iterator i = firstRates.begin(); i != firstRates.end(); i++)
		{
			for (std::vector<mpfr_class>::const_iterator j = secondRates.begin(); j != secondRates.end(); j++)
			{
				mpfr_class inverse = 1 / (*j - *i);
				scratchMemory[std::distance(firstRates.begin(), i)] *= inverse;
				scratchMemory[std::distance(secondRates.begin(), j) + firstRates.size()] *= -inverse;
			}
		}
		//Last little bit....
		for (std::vector<mpfr_class>::const_iterator i = secondRates.begin(); i != secondRates.end(); i++)
		{
			for (std::vector<mpfr_class>::const_iterator j = i+1; j != secondRates.end(); j++)
			{
				mpfr_class inverse = 1 / (*j - *i);
				scratchMemory[std::distance(secondRates.begin(), i) + firstRates.size()] *= inverse;
				scratchMemory[std::distance(secondRates.begin(), j) + firstRates.size()] *= -inverse;
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