#include "NetworkReliabilitySubObs.h"
#include "NetworkReliabilityObs.h"
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/math/distributions/binomial.hpp>
#include "graphAlgorithms.h"
namespace networkReliability
{
	NetworkReliabilitySubObs::NetworkReliabilitySubObs(Context const& context, boost::shared_array<EdgeState> state, int radius, int conditioningCount, conditioning_type conditioningProb)
		:context(context), state(state), radius(radius), conditioningCount(conditioningCount), fixedInop(0), conditioningProb(conditioningProb)
	{
		std::vector<int>& capacityVector = context.getCapacityVector();
		const std::size_t nEdges = boost::num_edges(context.getGraph());
		//Profiling has indicated that the push_back below costs in terms of calls to new[]. 
		int couldBeDeactivatedCounter = 0;
		for(int i = 0; i < nEdges; i++)
		{
			if(state[i] & (UNFIXED_INOP | UNFIXED_OP)) couldBeDeactivatedCounter++;
		}
		couldBeDeactivated.reserve(couldBeDeactivatedCounter);
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
	void NetworkReliabilitySubObs::getRadius1ReducedGraph(Context::internalGraph& outputGraph, int& minimumInoperative, std::vector<int>& edgeCounts, std::vector<int>& components, boost::detail::depth_first_visit_restricted_impl_helper<Context::internalGraph>::stackType& stack, std::vector<boost::default_color_type>& colorMap) const
	{
		int nComponents = countComponents(context, state.get(), components, stack, colorMap);
		edgeCounts.clear();
		edgeCounts.resize(nComponents * nComponents);

		//determine the rates between all the different super-vertices
		Context::internalGraph::edge_iterator current, end;
		const Context::internalGraph& graph = context.getGraph();
		boost::tie(current, end) = boost::edges(graph);
		for (; current != end; current++)
		{
			int edgeIndex = boost::get(boost::edge_index, graph, *current);
			//is it an edge between super-nodes?
			if (state[edgeIndex] & UNFIXED_MASK)
			{
				edgeCounts[components[current->m_source] + components[current->m_target] * nComponents]++;
			}
		}
		outputGraph = Context::internalGraph(nComponents);
		int edgeCounter = 0;
		for (int i = 0; i < nComponents; i++)
		{
			if (edgeCounts[i + i * nComponents] > 0)
			{
				boost::add_edge(i, i, edgeCounter++, outputGraph);
			}
			for (int j = i + 1; j < nComponents; j++)
			{
				if (edgeCounts[i + j * nComponents] + edgeCounts[j + i * nComponents] > 0)
				{
					boost::add_edge(i, j, edgeCounter++, outputGraph);
				}
			}
		}
		if (context.useMinCut())
		{
			minimumInoperative = minCut;
		}
		else
		{
			minimumInoperative = std::max(0, conditioningCount - fixedInop);
		}
	}
	void NetworkReliabilitySubObs::getRadius1ReducedGraphNoSelf(Context::internalGraph& outputGraph, int& minimumInoperative, std::vector<int>& edgeCounts, std::vector<int>& components, boost::detail::depth_first_visit_restricted_impl_helper<Context::internalGraph>::stackType& stack, std::vector<boost::default_color_type>& colorMap) const
	{
		int nComponents = countComponents(context, state.get(), components, stack, colorMap);
		edgeCounts.clear();
		edgeCounts.resize(nComponents * nComponents);

		//determine the rates between all the different super-vertices
		Context::internalGraph::edge_iterator current, end;
		const Context::internalGraph& graph = context.getGraph();
		boost::tie(current, end) = boost::edges(graph);
		for (; current != end; current++)
		{
			int edgeIndex = boost::get(boost::edge_index, graph, *current);
			//is it an edge between super-nodes?
			if (state[edgeIndex] & UNFIXED_MASK)
			{
				edgeCounts[components[current->m_source] + components[current->m_target] * nComponents]++;
			}
		}
		outputGraph = Context::internalGraph(nComponents);
		int edgeCounter = 0;
		for (int i = 0; i < nComponents; i++)
		{
			for (int j = i + 1; j < nComponents; j++)
			{
				if (edgeCounts[i + j * nComponents] + edgeCounts[j + i * nComponents] > 0)
				{
					boost::add_edge(i, j, edgeCounter++, outputGraph);
				}
			}
		}
		if (context.useMinCut())
		{
			minimumInoperative = minCut;
		}
		else
		{
			minimumInoperative = std::max(0, conditioningCount - fixedInop);
		}
	}
	const NetworkReliabilitySubObs::conditioning_type& NetworkReliabilitySubObs::getConditioningProb() const
	{
		return conditioningProb;
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
		std::size_t nDeactivated;
		if (context.useMinCut())
		{
			const ::TruncatedBinomialDistribution::TruncatedBinomialDistribution& distribution = context.getInopDistribution(std::max(minCut, conditioningCount - fixedInop), couldBeDeactivated.size(), couldBeDeactivated.size());
			nDeactivated = distribution(randomSource);
		}
		else
		{
			const ::TruncatedBinomialDistribution::TruncatedBinomialDistribution& distribution = context.getInopDistribution(std::max(0, conditioningCount - fixedInop), couldBeDeactivated.size(), couldBeDeactivated.size());
			nDeactivated = distribution(randomSource);
		}

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
		if (context.useMinCut())
		{
			conditioning_type newConditioningProb;
			int newConditioningCount = std::max(fixedInop + minCut, conditioningCount);
			if(fixedInop + minCut > conditioningCount && minCut > 0)
			{
				boost::math::binomial_distribution<> relevantBinomial((double)couldBeDeactivated.size(), context.getInoperationalProbabilityD());
				if(conditioningCount > fixedInop)
				{
					const ::TruncatedBinomialDistribution::TruncatedBinomialDistribution& relevantDistribution = context.getInopDistribution(conditioningCount - fixedInop, couldBeDeactivated.size(), couldBeDeactivated.size());
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
		else
		{
			return NetworkReliabilityObs(context, newState, conditioningCount, conditioningProb);
		}
	}
	NetworkReliabilitySubObs::conditioning_type NetworkReliabilitySubObs::getGeneratedObservationConditioningProb() const
	{
		if (context.useMinCut())
		{
			conditioning_type newConditioningProb;
			int newConditioningCount = std::max(fixedInop + minCut, conditioningCount);
			if (fixedInop + minCut > conditioningCount && minCut > 0)
			{
				boost::math::binomial_distribution<> relevantBinomial((double)couldBeDeactivated.size(), context.getInoperationalProbabilityD());
				if (conditioningCount > fixedInop)
				{
					const ::TruncatedBinomialDistribution::TruncatedBinomialDistribution& relevantDistribution = context.getInopDistribution(conditioningCount - fixedInop, couldBeDeactivated.size(), couldBeDeactivated.size());
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
			return newConditioningProb * conditioningProb;
		}
		return conditioningProb;
	}
	const EdgeState* NetworkReliabilitySubObs::getState() const
	{
		return state.get();
	}
	int NetworkReliabilitySubObs::getConditioningCount() const
	{
		return conditioningCount;
	}
	NetworkReliabilitySubObs::NetworkReliabilitySubObs(Context const& context)
		:context(context)
	{}
	NetworkReliabilitySubObs NetworkReliabilitySubObs::copyWithConditioningProb(const conditioning_type& conditioningProb) const
	{
		NetworkReliabilitySubObs copy(context);
		copy.conditioningProb = conditioningProb;
		copy.state = state;
		copy.radius = radius;
		copy.minCut = minCut;
		copy.couldBeDeactivated.insert(copy.couldBeDeactivated.begin(), couldBeDeactivated.begin(), couldBeDeactivated.end());
		copy.conditioningCount = conditioningCount;
		copy.fixedInop = fixedInop;
		return copy;
	}
	Context const& NetworkReliabilitySubObs::getContext() const
	{
		return context;
	}
}
