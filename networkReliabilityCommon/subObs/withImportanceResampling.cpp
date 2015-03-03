#include "subObs/withImportanceResampling.h"
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include "graphAlgorithms.h"
#include "seriesParallelReduction.hpp"
#include "subObs/getObservation.hpp"
#include "obs/getSubObservation.hpp"
#include <boost/range/algorithm/random_shuffle.hpp>
#include <boost/random/random_number_generator.hpp>
#include "obs/withSub.h"
namespace networkReliability
{
	namespace subObs
	{
		void withImportanceResampling::initialise()
		{
			std::vector<int>& capacityVector = context.capacityVector;
			const std::size_t nEdges = boost::num_edges(context.getGraph());
			//Profiling has indicated that the push_back below costs in terms of calls to new[]. 
			int unknownStateCounter = 0;
			for(std::size_t i = 0; i < nEdges; i++)
			{
				if(state[i] & (UNFIXED_INOP | UNFIXED_OP)) unknownStateCounter++;
			}
			unknownState.reserve(unknownStateCounter);
			for(std::size_t i = 0; i < nEdges; i++)
			{
				if(state[i] == FIXED_OP) 
				{
					capacityVector[2*i] = capacityVector[2*i + 1] = HIGH_CAPACITY;
				}
				else if(state[i] & (UNFIXED_INOP | UNFIXED_OP))
				{
					unknownState.push_back((int)i);
					capacityVector[2*i] = capacityVector[2*i + 1] = 1;
				}
				else 
				{
					capacityVector[2*i] = capacityVector[2*i + 1] = 0;
					fixedInop++;
				}
			}
			minCut = context.getMinCut(capacityVector);
			if (minCut < HIGH_CAPACITY)
			{
				mpfr_class newConditioningProb;
				generatedObservationConditioningCount = std::max(fixedInop + minCut, conditioningCount);
				if(fixedInop + minCut > conditioningCount && minCut > 0)
				{
					if(conditioningCount > fixedInop)
					{
						const ::TruncatedBinomialDistribution::TruncatedBinomialDistribution& relevantDistribution = context.getInopDistribution(conditioningCount - fixedInop, unknownState.size(), unknownState.size());
						newConditioningProb = 1 - relevantDistribution.getCumulativeProbability(minCut - conditioningCount + fixedInop - 1);
					}
					else
					{
						boost::math::binomial_distribution<mpfr_class> relevantBinomial((double)unknownState.size(), context.getInoperationalProbabilityD());
						newConditioningProb = boost::math::cdf(boost::math::complement(relevantBinomial, minCut - 1));
					}
				}
				else newConditioningProb = 1;
				generatedObservationConditioningProb = newConditioningProb*conditioningProb;
			}
			else
			{
				generatedObservationConditioningProb = 0;
			}
			std::vector<bool> unknownStateBools(nEdges, false);
			for(std::size_t edgeCounter = 0; edgeCounter < nEdges; edgeCounter++)
			{
				if(state[edgeCounter] & UNFIXED_MASK) unknownStateBools[edgeCounter] = true;
			}
			identifyMinCutEdges(context.edgeResidualCapacityVector, context.capacityVector, state.get(), context.getDirectedGraph(), context.colorVector, unknownStateBools, importanceSamplingEdges, context.getInterestVertices()[0], context.getInterestVertices()[1]);
		}
		int withImportanceResampling::getMinCut() const
		{
			return minCut;
		}
		const mpfr_class& withImportanceResampling::getGeneratedObservationConditioningProb() const
		{
			return generatedObservationConditioningProb;
		}
		const mpfr_class& withImportanceResampling::getResamplingProb() const
		{
			return resamplingProb;
		}
		const mpfr_class& withImportanceResampling::getConditioningProb() const
		{
			return conditioningProb;
		}
		withImportanceResampling::withImportanceResampling(withImportanceResampling&& other)
			: ::networkReliability::subObs::subObs(static_cast<::networkReliability::subObs::subObs&&>(other))
		{
			minCut = other.minCut;
			unknownState.swap(other.unknownState);
			conditioningCount = other.conditioningCount;
			fixedInop = other.fixedInop;
			conditioningProb = other.conditioningProb;
			generatedObservationConditioningCount = other.generatedObservationConditioningCount;
			generatedObservationConditioningProb = other.generatedObservationConditioningProb;
			nextSmallerRadius = other.nextSmallerRadius;
			boundaryEdges.swap(other.boundaryEdges);
			importanceSamplingEdges.swap(other.importanceSamplingEdges);
			resamplingProb = other.resamplingProb;
		}
		withImportanceResampling& withImportanceResampling::operator=(withImportanceResampling&& other)
		{
			*static_cast<::networkReliability::subObs::subObs*>(this) = static_cast<::networkReliability::subObs::subObs&&>(other);
			minCut = other.minCut;
			unknownState.swap(other.unknownState);
			conditioningCount = other.conditioningCount;
			fixedInop = other.fixedInop;
			conditioningProb = other.conditioningProb;
			generatedObservationConditioningCount = other.generatedObservationConditioningCount;
			generatedObservationConditioningProb = other.generatedObservationConditioningProb;
			nextSmallerRadius = other.nextSmallerRadius;
			boundaryEdges.swap(other.boundaryEdges);
			importanceSamplingEdges.swap(other.importanceSamplingEdges);
			resamplingProb = other.resamplingProb;
			return *this;
		}
		withImportanceResampling::withImportanceResampling(Context const& context, boost::shared_array<EdgeState> state, double radius, ::networkReliability::subObs::withImportanceResamplingConstructorType& other)
			: ::networkReliability::subObs::subObs(context, state, radius), conditioningCount(other.conditioningCount), fixedInop(0), conditioningProb(other.conditioningProb), nextSmallerRadius(other.nextRadius), resamplingProb(other.resamplingProb)
		{
			boundaryEdges.swap(other.boundaryEdges);
			initialise();
		}
		void withImportanceResampling::getObservation(EdgeState* newState, boost::mt19937& randomSource, observationConstructorType& otherData) const
		{
			const std::size_t nEdges = context.getNEdges();
			memcpy(newState, state.get(), sizeof(EdgeState)*nEdges);
			if(radius == 0)
			{
				otherData.conditioningCount = conditioningCount;
				otherData.conditioningProb = 1;
				otherData.resamplingProb = resamplingProb;
				return;
			}
			for(std::size_t i = 0; i < unknownState.size(); i++)
			{
				newState[unknownState[i]] = UNFIXED_OP;
			}

			int minAdditionalDeactivated = std::max(minCut, conditioningCount - fixedInop);
			const ::TruncatedBinomialDistribution::TruncatedBinomialDistribution& distribution = context.getInopDistribution(minAdditionalDeactivated, unknownState.size(), unknownState.size());
			std::size_t nDeactivated = distribution(randomSource);

			int maximumNumberOfImportanceResamplingEdges = std::min(nDeactivated, importanceSamplingEdges.size());
			int minimumNumberOfImportanceResamplingEdges = std::max((int)0, (int)(nDeactivated - (unknownState.size() - boundaryEdges.size()) - (boundaryEdges.size() - importanceSamplingEdges.size())));
			boost::random::uniform_int_distribution<> resamplingCountDistribution(minimumNumberOfImportanceResamplingEdges, maximumNumberOfImportanceResamplingEdges);
			int numberResamplingEdges = resamplingCountDistribution(randomSource);
			//grab some scratch memory
			std::vector<int>& permuted = context.distanceVector;
			permuted.clear();
			permuted.insert(permuted.begin(), importanceSamplingEdges.begin(), importanceSamplingEdges.end());

			boost::random_number_generator<boost::mt19937> generator(randomSource);
			boost::random_shuffle(permuted, generator);
			for(int i = 0; i < numberResamplingEdges; i++)
			{
				newState[permuted[i]] = UNFIXED_INOP;
			}
			for(int i = numberResamplingEdges; i < (int)permuted.size(); i++)
			{
				newState[permuted[i]] = FIXED_OP;
			}

			int numberNonResamplingEdges = nDeactivated - numberResamplingEdges;
			permuted.clear();
			permuted.insert(permuted.begin(), unknownState.begin(), unknownState.end());
			boost::random_shuffle(permuted, generator);
			for(int i = 0; i < numberNonResamplingEdges; i++)
			{
				if(newState[*permuted.rbegin()] == UNFIXED_INOP)
				{
					i--;
				}
				else if(newState[*permuted.rbegin()] == FIXED_OP)
				{
					newState[*permuted.rbegin()] = UNFIXED_OP;
					i--;
				}
				else newState[*permuted.rbegin()] = UNFIXED_INOP;
				permuted.pop_back();
			}
			for(std::size_t i = 0; i < permuted.size(); i++)
			{
				if(newState[permuted[i]] == FIXED_OP)
				{
					newState[permuted[i]] = UNFIXED_OP;
				}
			}
			boost::shared_array<EdgeState> followingSubObsState(new EdgeState[nEdges]);
			::networkReliability::withSub::getSubObservation(nextSmallerRadius, followingSubObsState.get(), context, newState);
			int followingSubObsFixedInop = 0;
			int followingSubObsUnknownState = 0;
			int importanceResamplingInop = 0, importanceResamplingUnknown = 0;
			for(std::vector<int>::const_iterator i = importanceSamplingEdges.begin(); i != importanceSamplingEdges.end(); i++)
			{
				if(followingSubObsState[*i] == FIXED_INOP) importanceResamplingInop++;
				else if(followingSubObsState[*i] & UNFIXED_MASK) importanceResamplingUnknown++;
			}
			for(std::size_t i = 0; i < nEdges; i++)
			{
				if(followingSubObsState[i] == FIXED_INOP)
				{
					followingSubObsFixedInop++;
				}
				else if(followingSubObsState[i] & UNFIXED_MASK)
				{
					followingSubObsUnknownState++;
				}
			}
			int boundaryFixedInop = followingSubObsFixedInop - fixedInop;
			int boundaryFixedOp = (unknownState.size() - followingSubObsUnknownState) - boundaryFixedInop;
			mpfr_class operationalProbability = context.getOperationalProbability();
			mpfr_class inoperationalProbability = 1 - operationalProbability;
			mpfr_class usualDensity = boost::multiprecision::pow(operationalProbability, boundaryFixedOp) * boost::multiprecision::pow(inoperationalProbability, boundaryFixedInop);
			//Finish calculating the usualDensity.
			const TruncatedBinomialDistribution::TruncatedBinomialDistribution& dist = context.getInopDistribution(0, unknownState.size(), unknownState.size());
			/* The unoptimised code*/
			/*for(int i = followingSubObsFixedInop; i < minAdditionalDeactivated+fixedInop; i++)
			{
				mpfr_class lower = 0;
				if(i - fixedInop > 0) lower = dist.getCumulativeProbability(i - fixedInop - 1);
				mpfr_class upper = dist.getCumulativeProbability(i - fixedInop);
				usualDensity -= (upper - lower) * boost::math::binomial_coefficient<mpfr_class>(followingSubObsUnknownState, i - followingSubObsFixedInop) / boost::math::binomial_coefficient<mpfr_class>(unknownState.size(), i-fixedInop);
			}
			mpfr_class importanceDensity = 0;
			for(int i = std::max(minAdditionalDeactivated, boundaryFixedInop); i < followingSubObsUnknownState + boundaryFixedInop + 1; i++)
			{
				mpfr_class currentPart = 0;
				int upperImportanceEdges = std::min(importanceResamplingInop + importanceResamplingUnknown, i - (boundaryFixedInop - importanceResamplingInop));
				int lowerImportanceEdges = std::max(i - (followingSubObsUnknownState - importanceResamplingUnknown) - (boundaryFixedInop - importanceResamplingInop), importanceResamplingInop);
				float weight = 1.0 / (std::min((int)importanceSamplingEdges.size(), i) - std::max(0, i - (int)(unknownState.size() - importanceSamplingEdges.size())) + 1);
				for(int j = lowerImportanceEdges; j < upperImportanceEdges + 1; j++)
				{
					currentPart += (boost::math::binomial_coefficient<mpfr_class>(importanceResamplingUnknown, j - importanceResamplingInop) / boost::math::binomial_coefficient<mpfr_class>(importanceSamplingEdges.size(), j)) * (boost::math::binomial_coefficient<mpfr_class>(followingSubObsUnknownState - importanceResamplingUnknown, i - j - (boundaryFixedInop - importanceResamplingInop))/boost::math::binomial_coefficient<mpfr_class>(unknownState.size()-importanceSamplingEdges.size(), i - j));
				}
				mpfr_class lower = 0;
				if(i > 0) lower = dist.getCumulativeProbability(i - 1);
				importanceDensity += weight * currentPart * (dist.getCumulativeProbability(i) - lower);
			}*/
			/* Optimised code */
			mpfr_class lower = 0;
			if(followingSubObsFixedInop < minAdditionalDeactivated+fixedInop && nextSmallerRadius > 0)
			{
				int i = followingSubObsFixedInop;
				mpfr_class binomRatio = boost::math::binomial_coefficient<mpfr_class>(followingSubObsUnknownState, i - followingSubObsFixedInop) / boost::math::binomial_coefficient<mpfr_class>(unknownState.size(), i-fixedInop);
				//double binomRatio = boost::math::binomial_coefficient<double>(followingSubObsUnknownState, i - followingSubObsFixedInop) / boost::math::binomial_coefficient<double>(unknownState.size(), i - fixedInop);
				if(followingSubObsFixedInop - fixedInop > 0) lower = dist.getCumulativeProbability(i - fixedInop - 1);
				mpfr_class upper = dist.getCumulativeProbability(i - fixedInop);
				usualDensity -= (upper - lower)*binomRatio;
				i++;
				for(; i < minAdditionalDeactivated+fixedInop; i++)
				{
					lower = 0;
					if(i - fixedInop > 0) lower = dist.getCumulativeProbability(i - fixedInop - 1);
					upper = dist.getCumulativeProbability(i - fixedInop);
					//binomRatio = boost::math::binomial_coefficient<mpfr_class>(followingSubObsUnknownState, i - followingSubObsFixedInop) / boost::math::binomial_coefficient<mpfr_class>(unknownState.size(), i-fixedInop);
					binomRatio *= (i - fixedInop) * (followingSubObsUnknownState - (i - followingSubObsFixedInop - 1))/(mpfr_class)((unknownState.size() - (i - fixedInop - 1))*(i - followingSubObsFixedInop));
					usualDensity -= (upper - lower) * binomRatio;
				}
			}
			mpfr_class importanceDensity = 0;
			if(std::max(minAdditionalDeactivated, boundaryFixedInop) < followingSubObsUnknownState + boundaryFixedInop + 1)
			{
				int i = std::max(minAdditionalDeactivated, boundaryFixedInop);
				mpfr_class currentPart = 0;
				int upperImportanceEdges = std::min(importanceResamplingInop + importanceResamplingUnknown, i - (boundaryFixedInop - importanceResamplingInop));
				int lowerImportanceEdges = std::max(i - (followingSubObsUnknownState - importanceResamplingUnknown) - (boundaryFixedInop - importanceResamplingInop), importanceResamplingInop);
				double weight = 1.0 / (std::min((int)importanceSamplingEdges.size(), i) - std::max(0, i - (int)(unknownState.size() - importanceSamplingEdges.size())) + 1);
				std::vector<mpfr_class> binomialRatios(upperImportanceEdges + 1, std::numeric_limits<double>::quiet_NaN());
				int j = lowerImportanceEdges;
				binomialRatios[j] = (boost::math::binomial_coefficient<mpfr_class>(importanceResamplingUnknown, j - importanceResamplingInop) / boost::math::binomial_coefficient<mpfr_class>(importanceSamplingEdges.size(), j))*(boost::math::binomial_coefficient<mpfr_class>(followingSubObsUnknownState - importanceResamplingUnknown, i - j - (boundaryFixedInop - importanceResamplingInop))/boost::math::binomial_coefficient<mpfr_class>(unknownState.size()-importanceSamplingEdges.size(), i - j));
				currentPart += binomialRatios[j];
				j++;
				for(; j < upperImportanceEdges + 1; j++)
				{
					binomialRatios[j] = binomialRatios[j-1]* (j*(importanceResamplingUnknown - ((j-1) - importanceResamplingInop))/(mpfr_class)((importanceSamplingEdges.size()-(j-1))*(j - importanceResamplingInop)))*((i - (j-1)- (boundaryFixedInop - importanceResamplingInop))*(unknownState.size()-importanceSamplingEdges.size()-(i - j))/(mpfr_class)((i-(j-1))*(followingSubObsUnknownState - importanceResamplingUnknown - (i - j - (boundaryFixedInop - importanceResamplingInop)))));
					currentPart += binomialRatios[j];
				}
				lower = 0;
				if(i > 0) lower = dist.getCumulativeProbability(i - 1);
				importanceDensity += weight * currentPart * (dist.getCumulativeProbability(i) - lower);
				i++;
				for(; i < followingSubObsUnknownState + boundaryFixedInop + 1; i++)
				{
					currentPart = 0;
					upperImportanceEdges = std::min(importanceResamplingInop + importanceResamplingUnknown, i - (boundaryFixedInop - importanceResamplingInop));
					lowerImportanceEdges = std::max(i - (followingSubObsUnknownState - importanceResamplingUnknown) - (boundaryFixedInop - importanceResamplingInop), importanceResamplingInop);
					weight = 1.0 / (std::min((int)importanceSamplingEdges.size(), i) - std::max(0, i - (int)(unknownState.size() - importanceSamplingEdges.size())) + 1);
					for(int j = lowerImportanceEdges; j < upperImportanceEdges + 1; j++)
					{
						if(j < (int)binomialRatios.size() && binomialRatios[j] == binomialRatios[j]) binomialRatios[j] *= (i-j)*(followingSubObsUnknownState - importanceResamplingUnknown - (i - j - (boundaryFixedInop - importanceResamplingInop) - 1))/(mpfr_class)((i - j - (boundaryFixedInop - importanceResamplingInop))*(unknownState.size()-importanceSamplingEdges.size() - (i - j - 1)));
						else
						{
							if(j >= (int)binomialRatios.size()) binomialRatios.resize(j+1);
							binomialRatios[j] = (boost::math::binomial_coefficient<mpfr_class>(importanceResamplingUnknown, j - importanceResamplingInop) / boost::math::binomial_coefficient<mpfr_class>(importanceSamplingEdges.size(), j)) * (boost::math::binomial_coefficient<mpfr_class>(followingSubObsUnknownState - importanceResamplingUnknown, i - j - (boundaryFixedInop - importanceResamplingInop))/boost::math::binomial_coefficient<mpfr_class>(unknownState.size()-importanceSamplingEdges.size(), i - j));
						}
						currentPart += binomialRatios[j];
					}
					lower = 0;
					if(i > 0) lower = dist.getCumulativeProbability(i - 1);
					importanceDensity += weight * currentPart * (dist.getCumulativeProbability(i) - lower);
				}
			}
			otherData.conditioningCount = generatedObservationConditioningCount;
			otherData.conditioningProb = generatedObservationConditioningProb;
			otherData.resamplingProb = resamplingProb * usualDensity/importanceDensity;
		}
		int withImportanceResampling::getConditioningCount() const
		{
			return conditioningCount;
		}
		withImportanceResampling::withImportanceResampling(Context const& context, boost::shared_array<EdgeState> state, double radius)
			: ::networkReliability::subObs::subObs(context, state, radius)
		{}
		withImportanceResampling withImportanceResampling::copyWithGeneratedObservationConditioningProb(const mpfr_class& newGeneratedObservationConditioningProb) const
		{
			withImportanceResampling copy(context, state, radius);
			copy.conditioningProb = conditioningProb;
			copy.generatedObservationConditioningProb = newGeneratedObservationConditioningProb;
			copy.minCut = minCut;
			copy.unknownState.insert(copy.unknownState.begin(), unknownState.begin(), unknownState.end());
			copy.conditioningCount = conditioningCount;
			copy.generatedObservationConditioningCount = generatedObservationConditioningCount;
			copy.fixedInop = fixedInop;
			copy.nextSmallerRadius = nextSmallerRadius;
			copy.boundaryEdges = boundaryEdges;
			copy.importanceSamplingEdges = importanceSamplingEdges;
			copy.resamplingProb = resamplingProb;
			return copy;
		}
		int withImportanceResampling::getFixedInopCount() const
		{
			return fixedInop;
		}
		const std::vector<int>& withImportanceResampling::getUnknownStateEdges() const
		{
			return unknownState;
		}
	}
}
