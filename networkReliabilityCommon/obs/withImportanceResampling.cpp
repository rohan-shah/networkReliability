#include "obs/withImportanceResampling.h"
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/math/distributions/binomial.hpp>
namespace networkReliability
{
	namespace obs
	{
		withImportanceResampling::withImportanceResampling(Context const& context, boost::mt19937& randomSource)
			: ::networkReliability::withSub(context, randomSource), conditioningCount(0), conditioningProb(1) 
		{}
		withImportanceResampling withImportanceResampling::constructConditional(Context const& context, boost::mt19937& randomSource)
		{
			const Context::internalGraph& graph = context.getGraph();
			const std::size_t nEdges = boost::num_edges(graph);
			const std::size_t minCutEdges = context.getMinCutEdges();
			boost::shared_array<EdgeState> state(new EdgeState[nEdges]);
			::networkReliability::NetworkReliabilityObs::constructConditional(context, randomSource, state.get(), false);

			const ::TruncatedBinomialDistribution::TruncatedBinomialDistribution& originalDist = context.getOpDistribution(0, nEdges, nEdges);
			mpfr_class conditioningProb = originalDist.getCumulativeProbability((int)(nEdges - minCutEdges));
			return withImportanceResampling(context, state, (int)minCutEdges, conditioningProb);
		}
		withImportanceResampling::withImportanceResampling(Context const& context, boost::shared_array<EdgeState> state, int conditioningCount, mpfr_class conditioningProb)
			:  ::networkReliability::withSub(context, state), conditioningCount(conditioningCount), conditioningProb(conditioningProb)
		{}
		withImportanceResampling& withImportanceResampling::operator=(withImportanceResampling&& other)
		{
			this->::networkReliability::withSub::operator=(static_cast<::networkReliability::withSub&&>(other));
			conditioningCount = other.conditioningCount;
			conditioningProb = other.conditioningProb;
			return *this;
		}
		withImportanceResampling::withImportanceResampling(withImportanceResampling&& other)
			: ::networkReliability::withSub(static_cast<withImportanceResampling&&>(other)), conditioningCount(other.conditioningCount), conditioningProb(other.conditioningProb)
		{}
		mpfr_class withImportanceResampling::getConditioningProb() const
		{
			return conditioningProb;
		}
		int withImportanceResampling::getConditioningCount() const
		{
			return conditioningCount;
		}
	}
}
