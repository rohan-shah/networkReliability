#include "subObs/withImportanceResampling.h"
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/math/distributions/binomial.hpp>
#include "graphAlgorithms.h"
#include "seriesParallelReduction.hpp"
#include "subObs/getObservation.hpp"
#include "obs/getSubObservation.hpp"
namespace networkReliability
{
	namespace subObs
	{
		withImportanceResampling::withImportanceResampling(Context const& context, boost::shared_array<EdgeState> state, double radius, int conditioningCount, mpfr_class conditioningProb)
			: ::networkReliability::subObs::subObs(context, state, radius), conditioningCount(conditioningCount), fixedInop(0), conditioningProb(conditioningProb)
		{
			std::vector<int>& capacityVector = context.capacityVector;
			const std::size_t nEdges = boost::num_edges(context.getGraph());
			//Profiling has indicated that the push_back below costs in terms of calls to new[]. 
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
					couldBeDeactivated.push_back((int)i);
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
		int withImportanceResampling::getMinCut() const
		{
			return minCut;
		}
		const mpfr_class& withImportanceResampling::getGeneratedObservationConditioningProb() const
		{
			return generatedObservationConditioningProb;
		}
		const mpfr_class& withImportanceResampling::getConditioningProb() const
		{
			return conditioningProb;
		}
		withImportanceResampling::withImportanceResampling(withImportanceResampling&& other)
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
		withImportanceResampling& withImportanceResampling::operator=(withImportanceResampling&& other)
		{
			*static_cast<::networkReliability::subObs::subObs*>(this) = static_cast<::networkReliability::subObs::subObs&&>(other);
			minCut = other.minCut;
			couldBeDeactivated.swap(other.couldBeDeactivated);
			conditioningCount = other.conditioningCount;
			fixedInop = other.fixedInop;
			conditioningProb = other.conditioningProb;
			generatedObservationConditioningCount = other.generatedObservationConditioningCount;
			generatedObservationConditioningProb = other.generatedObservationConditioningProb;
			return *this;
		}
		void withImportanceResampling::getObservation(EdgeState* newState, boost::mt19937& randomSource, observationConstructorType& otherData) const
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
			copy.couldBeDeactivated.insert(copy.couldBeDeactivated.begin(), couldBeDeactivated.begin(), couldBeDeactivated.end());
			copy.conditioningCount = conditioningCount;
			copy.generatedObservationConditioningCount = generatedObservationConditioningCount;
			copy.fixedInop = fixedInop;
			return copy;
		}
		int withImportanceResampling::getFixedInopCount() const
		{
			return fixedInop;
		}
		const std::vector<int>& withImportanceResampling::getPotentiallyDeactivated() const
		{
			return couldBeDeactivated;
		}
	}
}
