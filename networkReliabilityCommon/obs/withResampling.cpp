#include "obs/withResampling.h"
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/math/distributions/binomial.hpp>
namespace networkReliability
{
	namespace obs
	{
		withResampling::withResampling(context const& contextObj, boost::mt19937& randomSource)
			: ::networkReliability::withSub(contextObj, randomSource), conditioningCount(0), conditioningProb(1)
		{
		}
		withResampling withResampling::constructConditional(context const& contextObj, boost::mt19937& randomSource)
		{
			const context::internalGraph& graph = contextObj.getGraph();
			const std::size_t nEdges = boost::num_edges(graph);
			const std::size_t minCutEdges = contextObj.getMinCutEdges();
			boost::shared_array<edgeState> state(new edgeState[nEdges]);
			::networkReliability::NetworkReliabilityObs::constructConditional(contextObj, randomSource, state.get(), false);

			const ::TruncatedBinomialDistribution::TruncatedBinomialDistribution& originalDist = contextObj.getOpDistribution(0, nEdges, nEdges);
			mpfr_class conditioningProb = originalDist.getCumulativeProbability((int)nEdges - (int)minCutEdges);
			return withResampling(contextObj, state, (int)minCutEdges, conditioningProb);
		}
		withResampling::withResampling(context const& contextObj, boost::shared_array<edgeState> state, ::networkReliability::obs::withResamplingConstructorType& other)
			: ::networkReliability::withSub(contextObj, state), conditioningCount(other.conditioningCount), conditioningProb(other.conditioningProb)
		{
		}
		withResampling::withResampling(context const& contextObj, boost::shared_array<edgeState> state, int conditioningCount, mpfr_class conditioningProb)
			: ::networkReliability::withSub(contextObj, state), conditioningCount(conditioningCount), conditioningProb(conditioningProb)
		{}
		withResampling& withResampling::operator=(withResampling&& other)
		{
			this->withSub::operator=(static_cast<::networkReliability::withSub&&>(other));
			conditioningCount = other.conditioningCount;
			conditioningProb = other.conditioningProb;
			return *this;
		}
		withResampling::withResampling(withResampling&& other)
			: ::networkReliability::withSub(static_cast<::networkReliability::withSub&&>(other)), conditioningCount(other.conditioningCount), conditioningProb(other.conditioningProb)
		{}
		const mpfr_class& withResampling::getConditioningProb() const
		{
			return conditioningProb;
		}
		void withResampling::getSubObservation(edgeState* newStates, double radius, subObservationConstructorType& otherData) const
		{
			::networkReliability::withSub::getSubObservation(radius, newStates);
			otherData.conditioningCount = conditioningCount;
			otherData.conditioningProb = conditioningProb;
		}
		int withResampling::getConditioningCount() const
		{
			return conditioningCount;
		}
	}
}
