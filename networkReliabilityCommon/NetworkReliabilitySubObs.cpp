#include "NetworkReliabilitySubObs.h"
#include "NetworkReliabilityObs.h"
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/math/distributions/binomial.hpp>
namespace networkReliability
{
	NetworkReliabilitySubObs::NetworkReliabilitySubObs(Context const& context, boost::shared_array<EdgeState> state, int radius, int conditioningCount, conditioning_type conditioningProb)
		:context(context), state(state), radius(radius), conditioningCount(conditioningCount), fixedInop(0), conditioningProb(conditioningProb)
	{
		std::vector<int>& capacityVector = context.getCapacityVector();
		const std::size_t nEdges = boost::num_edges(context.getGraph());
		for(int i = 0; i < nEdges; i++)
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
		const Context::internalDirectedGraph& directedGraph = context.getDirectedGraph();
		minCut = context.getMinCut(capacityVector);
	}
	int NetworkReliabilitySubObs::getMinCut() const
	{
		return minCut;
	}
	NetworkReliabilitySubObs::NetworkReliabilitySubObs(NetworkReliabilitySubObs&& other)
		:context(other.context)
	{
		state.swap(other.state);
		radius = other.radius;
		minCut = other.minCut;
		couldBeDeactivated.swap(other.couldBeDeactivated);
		conditioningCount = other.conditioningCount;
		fixedInop = other.fixedInop;
		conditioningProb = other.conditioningProb;
	}
	NetworkReliabilitySubObs& NetworkReliabilitySubObs::operator=(NetworkReliabilitySubObs&& other)
	{
		state.swap(other.state);
		radius = other.radius;
		minCut = other.minCut;
		couldBeDeactivated.swap(other.couldBeDeactivated);
		conditioningCount = other.conditioningCount;
		fixedInop = other.fixedInop;
		conditioningProb = other.conditioningProb;
		return *this;
	}
	NetworkReliabilityObs NetworkReliabilitySubObs::getObservation(boost::mt19937& randomSource) const
	{
		if(radius == 0)
		{
			return NetworkReliabilityObs(context, state, conditioningCount, 1.0);
		}
		const std::size_t nEdges = context.getNEdges();
		//The minimum and maximum number of deactivated edges (before we do the conditioning step which says we need fixedInop + minCut to reach the end
		const ::TruncatedBinomialDistribution::TruncatedBinomialDistribution& distribution = context.getDistribution(std::max(minCut, conditioningCount - fixedInop), couldBeDeactivated.size(), couldBeDeactivated.size());
		std::size_t nDeactivated = distribution(randomSource);

		boost::shared_array<EdgeState> newState(new EdgeState[nEdges]);
		memcpy(newState.get(), state.get(), sizeof(EdgeState)*nEdges);

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
		conditioning_type newConditioningProb;
		int newConditioningCount = std::max(fixedInop + minCut, conditioningCount);
		if(fixedInop + minCut > conditioningCount && minCut > 0)
		{
			boost::math::binomial_distribution<> relevantBinomial((double)couldBeDeactivated.size(), context.getInoperationalProbabilityD());
			if(conditioningCount > fixedInop)
			{
				const ::TruncatedBinomialDistribution::TruncatedBinomialDistribution& relevantDistribution = context.getDistribution(conditioningCount - fixedInop, couldBeDeactivated.size(), couldBeDeactivated.size());
				const conditioning_type* cdf = relevantDistribution.getCumulativeProbabilities();
				newConditioningProb = 1 - cdf[minCut - conditioningCount + fixedInop - 1];

				/*double tmp1 = (1 - boost::math::cdf(relevantBinomial, minCut - 1)) / (1 - boost::math::cdf(relevantBinomial, conditioningCount - fixedInop - 1));
				double tmp2 = 1 - cdf[minCut - conditioningCount + fixedInop - 1].get_d();
				assert(abs(tmp1 - tmp2) < 1e-6);
				newConditioningProb = conditioning_type(1 - boost::math::cdf(relevantBinomial, minCut - 1)) / conditioning_type(1 - boost::math::cdf(relevantBinomial, conditioningCount - fixedInop - 1));*/
			}
			else
			{
				newConditioningProb = (1 - boost::math::cdf(relevantBinomial, minCut - 1));
			}
		}
		else newConditioningProb = 1;
		return NetworkReliabilityObs(context, newState, newConditioningCount, newConditioningProb*conditioningProb);
	}
	const EdgeState* NetworkReliabilitySubObs::getState() const
	{
		return state.get();
	}
	int NetworkReliabilitySubObs::getConditioningCount() const
	{
		return conditioningCount;
	}
}