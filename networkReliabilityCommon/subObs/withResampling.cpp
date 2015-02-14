#include "subObs/withResampling.h"
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/math/distributions/binomial.hpp>
#include "graphAlgorithms.h"
namespace networkReliability
{
	namespace subObs
	{
		withResampling::withResampling(Context const& context, boost::shared_array<EdgeState> state, double radius, int conditioningCount, mpfr_class conditioningProb)
		: ::networkReliability::subObs::subObs(context, state, radius), conditioningCount(conditioningCount), fixedInop(0), conditioningProb(conditioningProb)
		{
			initialise();
		}
		void withResampling::initialise()
		{
			const std::size_t nEdges = context.getNEdges();
			std::vector<int>& capacityVector = context.capacityVector;
			int couldBeDeactivatedCounter = 0;
			for(std::size_t i = 0; i < nEdges; i++)
			{
				if(state[i] & (UNFIXED_INOP | UNFIXED_OP)) couldBeDeactivatedCounter++;
			}
			couldBeDeactivated.reserve(couldBeDeactivatedCounter);
			for(std::size_t i = 0; i < nEdges; i++)
			{
				if(state[i] == FIXED_OP) 
				{
					capacityVector[2*i] = capacityVector[2*i + 1] = HIGH_CAPACITY;
				}
				else if(state[i] & (UNFIXED_INOP | UNFIXED_OP))
				{
					couldBeDeactivated.push_back(i);
					capacityVector[2*i] = capacityVector[2*i + 1] = 1;
				}
				else 
				{
					capacityVector[2*i] = capacityVector[2*i + 1] = 0;
					fixedInop++;
				}
			}

			minCut = context.getMinCut(capacityVector);
			if(minCut >= HIGH_CAPACITY)
			{}
			else
			{
				mpfr_class newConditioningProb;
				generatedObservationConditioningCount = std::max(fixedInop + minCut, conditioningCount);
				if(fixedInop + minCut > conditioningCount && minCut > 0)
				{
					if(conditioningCount > fixedInop)
					{
						const ::TruncatedBinomialDistribution::TruncatedBinomialDistribution& relevantDistribution = context.getInopDistribution(conditioningCount - fixedInop, couldBeDeactivated.size(), couldBeDeactivated.size());
						newConditioningProb = 1 - relevantDistribution.getCumulativeProbability(minCut - conditioningCount + fixedInop - 1);
					}
					else
					{
						boost::math::binomial_distribution<mpfr_class> relevantBinomial((double)couldBeDeactivated.size(), context.getInoperationalProbabilityD());
						newConditioningProb = boost::math::cdf(boost::math::complement(relevantBinomial, minCut - 1));
					}
				}
				else newConditioningProb = 1;
				generatedObservationConditioningProb = newConditioningProb*conditioningProb;
			}
		}
		int withResampling::getMinCut() const
		{
			return minCut;
		}
		const mpfr_class& withResampling::getGeneratedObservationConditioningProb() const
		{
			return generatedObservationConditioningProb;
		}
		const mpfr_class& withResampling::getConditioningProb() const
		{
			return conditioningProb;
		}
		withResampling::withResampling(withResampling&& other)
			: ::networkReliability::subObs::subObs(static_cast<::networkReliability::subObs::subObs&&>(other))
		{
			minCut = other.minCut;
			couldBeDeactivated.swap(other.couldBeDeactivated);
			conditioningCount = other.conditioningCount;
			fixedInop = other.fixedInop;
			conditioningProb = other.conditioningProb;
			generatedObservationConditioningCount = other.generatedObservationConditioningCount;
			generatedObservationConditioningProb = other.generatedObservationConditioningProb;
		}
		withResampling& withResampling::operator=(withResampling&& other)
		{
			this->::networkReliability::subObs::subObs::operator=(static_cast<::networkReliability::subObs::subObs&&>(other));
			minCut = other.minCut;
			couldBeDeactivated.swap(other.couldBeDeactivated);
			conditioningCount = other.conditioningCount;
			fixedInop = other.fixedInop;
			conditioningProb = other.conditioningProb;
			generatedObservationConditioningCount = other.generatedObservationConditioningCount;
			generatedObservationConditioningProb = other.generatedObservationConditioningProb;
			return *this;
		}
		void withResampling::getObservation(EdgeState* newState, boost::mt19937& randomSource, observationConstructorType& otherData) const
		{
			const std::size_t nEdges = context.getNEdges();
			memcpy(newState, state.get(), sizeof(EdgeState)*nEdges);
			if(radius == 0)
			{
				otherData.conditioningCount = conditioningCount;
				otherData.conditioningProb = 1;
				return;
			}
			const ::TruncatedBinomialDistribution::TruncatedBinomialDistribution& distribution = context.getInopDistribution(std::max(minCut, conditioningCount - fixedInop), couldBeDeactivated.size(), couldBeDeactivated.size());
			std::size_t nDeactivated = distribution(randomSource);

			for(std::size_t i = 0; i < couldBeDeactivated.size(); i++)
			{
				newState[couldBeDeactivated[i]] = UNFIXED_OP;
			}
			std::vector<int> indices(boost::counting_iterator<int>(0), boost::counting_iterator<int>((int)couldBeDeactivated.size()));
			for(std::size_t i = 0; i < nDeactivated; i++)
			{
				boost::random::uniform_int_distribution<std::size_t> dist(0, indices.size()-1);
				std::size_t generated = dist(randomSource);
				newState[couldBeDeactivated[indices[generated]]] = UNFIXED_INOP;
				std::swap(indices[generated], *indices.rbegin());
				indices.pop_back();
			}
			otherData.conditioningCount = generatedObservationConditioningCount;
			otherData.conditioningProb = generatedObservationConditioningProb;
		}
		int withResampling::getConditioningCount() const
		{
			return conditioningCount;
		}
		withResampling::withResampling(Context const& context, boost::shared_array<EdgeState> state, double radius)
			: ::networkReliability::subObs::subObs(context, state, radius)
		{
			//This is only called in copyWithGeneratedObservationConditioningProb, so don't initialise anything else
		}
		withResampling withResampling::copyWithGeneratedObservationConditioningProb(const mpfr_class& newGeneratedObservationConditioningProb) const
		{
			withResampling copy(context, state, radius);
			copy.conditioningProb = conditioningProb;
			copy.generatedObservationConditioningProb = newGeneratedObservationConditioningProb;
			copy.minCut = minCut;
			copy.couldBeDeactivated.insert(copy.couldBeDeactivated.begin(), couldBeDeactivated.begin(), couldBeDeactivated.end());
			copy.conditioningCount = conditioningCount;
			copy.generatedObservationConditioningCount = generatedObservationConditioningCount;
			copy.fixedInop = fixedInop;
			return copy;
		}
		int withResampling::getFixedInopCount() const
		{
			return fixedInop;
		}
		const std::vector<int>& withResampling::getPotentiallyDeactivated() const
		{
			return couldBeDeactivated;
		}
		withResampling::withResampling(Context const& context, boost::shared_array<EdgeState> state, double radius, ::networkReliability::subObs::withResamplingConstructorType& other)
			: ::networkReliability::subObs::subObs(context, state, radius), conditioningCount(other.conditioningCount), fixedInop(0), conditioningProb(other.conditioningProb)
		{
			initialise();
		}
	}
}
