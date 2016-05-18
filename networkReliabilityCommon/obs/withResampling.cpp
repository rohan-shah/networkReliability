#include "obs/withResampling.h"
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/math/distributions/binomial.hpp>
namespace networkReliability
{
	namespace obs
	{
		withResampling::withResampling(Context const& context, boost::mt19937& randomSource)
			: ::networkReliability::withSub(context, randomSource), conditioningCount(0), conditioningProb(1)
		{
		}
		withResampling withResampling::constructConditional(Context const& context, boost::mt19937& randomSource)
		{
			const Context::internalGraph& graph = context.getGraph();
			const std::size_t nEdges = boost::num_edges(graph);
			const std::size_t minCutEdges = context.getMinCutEdges();
			boost::shared_array<edgeState> state(new edgeState[nEdges]);
			::networkReliability::NetworkReliabilityObs::constructConditional(context, randomSource, state.get(), false);

			const ::TruncatedBinomialDistribution::TruncatedBinomialDistribution& originalDist = context.getOpDistribution(0, nEdges, nEdges);
			mpfr_class conditioningProb = originalDist.getCumulativeProbability((int)nEdges - (int)minCutEdges);
			return withResampling(context, state, (int)minCutEdges, conditioningProb);
		}
		withResampling::withResampling(Context const& context, boost::shared_array<edgeState> state, ::networkReliability::obs::withResamplingConstructorType& other)
			: ::networkReliability::withSub(context, state), conditioningCount(other.conditioningCount), conditioningProb(other.conditioningProb)
		{
		}
		withResampling::withResampling(Context const& context, boost::shared_array<edgeState> state, int conditioningCount, mpfr_class conditioningProb)
			: ::networkReliability::withSub(context, state), conditioningCount(conditioningCount), conditioningProb(conditioningProb)
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
