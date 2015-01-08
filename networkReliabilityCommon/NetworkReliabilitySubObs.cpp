#include "NetworkReliabilitySubObs.h"
#include "NetworkReliabilityObs.h"
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/math/distributions/binomial.hpp>
#include "graphAlgorithms.h"
#include "seriesParallelReduction.hpp"
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
		if(minCut >= HIGH_CAPACITY)
		{}
		else if (context.useMinCut())
		{
			conditioning_type newConditioningProb;
			generatedObservationConditioningCount = std::max(fixedInop + minCut, conditioningCount);
			if(fixedInop + minCut > conditioningCount && minCut > 0)
			{
				if(conditioningCount > fixedInop)
				{
					const ::TruncatedBinomialDistribution::TruncatedBinomialDistribution& relevantDistribution = context.getInopDistribution(conditioningCount - fixedInop, couldBeDeactivated.size(), couldBeDeactivated.size());
					const conditioning_type* cdf = relevantDistribution.getCumulativeProbabilities();
					newConditioningProb = 1 - cdf[minCut - conditioningCount + fixedInop - 1];
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
		else
		{
			generatedObservationConditioningProb = conditioningProb;
			generatedObservationConditioningCount = conditioningCount;
		}

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
	void NetworkReliabilitySubObs::getRadius1ReducedGraphNoSelfWithWeights(getRadius1ReducedGraphNoSelfWithWeightsInput& input) const
	{
		input.nUnreducedEdges = 0;
		int nComponents = countComponents(context, state.get(), input.components, input.stack, input.colorMap);
		input.edgeCounts.clear();
		input.edgeCounts.resize(nComponents * nComponents);

		//determine the rates between all the different super-vertices
		Context::internalGraph::edge_iterator currentBaseGraph, endBaseGraph;
		const Context::internalGraph& graph = context.getGraph();
		boost::tie(currentBaseGraph, endBaseGraph) = boost::edges(graph);
		for (; currentBaseGraph != endBaseGraph; currentBaseGraph++)
		{
			int edgeIndex = boost::get(boost::edge_index, graph, *currentBaseGraph);
			//is it an edge between super-nodes?
			if (state[edgeIndex] & UNFIXED_MASK)
			{
				input.edgeCounts[input.components[currentBaseGraph->m_source] + input.components[currentBaseGraph->m_target] * nComponents]++;
				input.nUnreducedEdges++;
			}
		}
		//Construct graph where every connected chunk of vertices is now a single vertex
		input.outputGraph = reducedGraphWithProbabilities(nComponents);
		//Put in operational / inoperational probabilities for every edge
		int edgeCounter = 0;
		const mpfr_class& opProbability = context.getOperationalProbability();
		mpfr_class inopProbability = 1 - opProbability;
		for (int i = 0; i < nComponents; i++)
		{
			boost::put(boost::vertex_name, input.outputGraph, i, i);
			for (int j = i + 1; j < nComponents; j++)
			{
				if (input.edgeCounts[i + j * nComponents] + input.edgeCounts[j + i * nComponents] > 0)
				{
					std::pair<reducedGraphWithProbabilities::edge_descriptor, bool> addedEdge = boost::add_edge(i, j, edgeCounter++, input.outputGraph);
					mpfr_class thisEdgeInopProbability = boost::multiprecision::pow(inopProbability, input.edgeCounts[i + j * nComponents] + input.edgeCounts[j + i * nComponents]);
					mpfr_class thisEdgeOpProbability = 1 - thisEdgeInopProbability;
					boost::put(boost::edge_op_probability, input.outputGraph, addedEdge.first, thisEdgeOpProbability);
					boost::put(boost::edge_inop_probability, input.outputGraph, addedEdge.first, thisEdgeInopProbability);
				}
			}
		}
		//Work out which vertices are the interest vertices in the reduced graph
		input.reducedInterestVertices.clear();
		for(std::vector<int>::const_iterator i = input.interestVertices.begin(); i != input.interestVertices.end(); i++)
		{
			input.reducedInterestVertices.push_back(input.components[*i]);
		}
		//Do series / parallel reduction
		seriesParallelReduction(input.outputGraph, input.reducedInterestVertices);
		//The reduced graph now has a different number of vertices to what it started with, but the names of the reduced vertices still correspond to connected components of the original graph. So we work backwards and get out the interest vertices in the reduced graph
		input.edgeCounts.clear();
		input.reducedInterestVertices.clear();
		input.reducedInterestVertices.resize(input.interestVertices.size());
		reducedGraphWithProbabilities::vertex_iterator current, end;
		boost::tie(current, end) = boost::vertices(input.outputGraph);
		for(; current != end; current++)
		{
			for(int i = 0; i < input.interestVertices.size(); i++)
			{
				if(boost::get(boost::vertex_name, input.outputGraph, *current) == input.components[input.interestVertices[i]])
				{
					input.reducedInterestVertices[i] = *current;
				}
			}
		}
		//Reset the edge_indices to be unique consecutive integers. 
		reducedGraphWithProbabilities::edge_iterator currentEdge, endEdge;
		boost::tie(currentEdge, endEdge) = boost::edges(input.outputGraph);
		int counter = 0;
		for(;currentEdge != endEdge; currentEdge++, counter++)
		{
			boost::put(boost::edge_index, input.outputGraph, *currentEdge, counter);
		}
	}
	const NetworkReliabilitySubObs::conditioning_type& NetworkReliabilitySubObs::getGeneratedObservationConditioningProb() const
	{
		return generatedObservationConditioningProb;
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
		generatedObservationConditioningCount = other.generatedObservationConditioningCount;
		generatedObservationConditioningProb = other.generatedObservationConditioningProb;
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
		generatedObservationConditioningCount = other.generatedObservationConditioningCount;
		generatedObservationConditioningProb = other.generatedObservationConditioningProb;
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
		return NetworkReliabilityObs(context, newState, generatedObservationConditioningCount, generatedObservationConditioningProb);
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
	NetworkReliabilitySubObs NetworkReliabilitySubObs::copyWithGeneratedObservationConditioningProb(const conditioning_type& newGeneratedObservationConditioningProb) const
	{
		NetworkReliabilitySubObs copy(context);
		copy.conditioningProb = conditioningProb;
		copy.generatedObservationConditioningProb = newGeneratedObservationConditioningProb;
		copy.state = state;
		copy.radius = radius;
		copy.minCut = minCut;
		copy.couldBeDeactivated.insert(copy.couldBeDeactivated.begin(), couldBeDeactivated.begin(), couldBeDeactivated.end());
		copy.conditioningCount = conditioningCount;
		copy.generatedObservationConditioningCount = generatedObservationConditioningCount;
		copy.fixedInop = fixedInop;
		return copy;
	}
	Context const& NetworkReliabilitySubObs::getContext() const
	{
		return context;
	}
	int NetworkReliabilitySubObs::getRadius() const
	{
		return radius;
	}
}
