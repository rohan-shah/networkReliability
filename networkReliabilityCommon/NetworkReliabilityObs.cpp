#include "NetworkReliabilityObs.h"
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/math/distributions/binomial.hpp>
#include "graphAlgorithms.h"
#include "seriesParallelReduction.hpp"
namespace networkReliability
{
	NetworkReliabilityObs::NetworkReliabilityObs(Context const& context, boost::archive::binary_iarchive& archive)
		:context(context)
	{
		archive >> *this;
	}
	NetworkReliabilityObs::NetworkReliabilityObs(Context const& context, boost::archive::text_iarchive& archive)
		:context(context)
	{
		archive >> *this;
	}
	NetworkReliabilityObs::NetworkReliabilityObs(Context const& context, boost::mt19937& randomSource)
		:context(context)
	{
		const Context::internalGraph& graph = context.getGraph();
		const std::size_t nEdges = boost::num_edges(graph);
		boost::shared_array<EdgeState> state(new EdgeState[nEdges]);

		boost::random::bernoulli_distribution<double> edgeDist(context.getOperationalProbability().convert_to<double>());

		for(std::size_t i = 0; i < nEdges; i++)
		{
			if(edgeDist(randomSource))
			{
				state[i] = UNFIXED_OP;
			}
			else state[i] = UNFIXED_INOP;
		}
		this->state = state;
	}
	void NetworkReliabilityObs::constructConditional(Context const& context, boost::mt19937& randomSource, EdgeState* state, bool fixed)
	{
		const Context::internalGraph& graph = context.getGraph();
		const std::size_t nEdges = boost::num_edges(graph);
		const std::size_t minCutEdges = context.getMinCutEdges();
		const ::TruncatedBinomialDistribution::TruncatedBinomialDistribution& dist = context.getInopDistribution(minCutEdges, nEdges, nEdges);

		const std::size_t nRemovedEdges = dist(randomSource);
		EdgeState inopState;
		if(fixed) 
		{
			std::fill(state, state + nEdges, FIXED_OP);
			inopState = FIXED_INOP;
		}
		else 
		{
			std::fill(state, state + nEdges, UNFIXED_OP);
			inopState = UNFIXED_INOP;
		}

		std::vector<int> indices(boost::counting_iterator<int>(0), boost::counting_iterator<int>((int)nEdges));
		for(std::size_t i = 0; i < nRemovedEdges; i++)
		{
			boost::random::uniform_int_distribution<int> removedEdgeIndexDistribution(0, (int)indices.size()-1);
			int index = removedEdgeIndexDistribution(randomSource);
			state[indices[index]] = inopState;
			std::swap(indices[index], *indices.rbegin());
			indices.pop_back();
		}
	}
	NetworkReliabilityObs::NetworkReliabilityObs(Context const& context, boost::shared_array<EdgeState> state)
		:context(context), state(state)
	{}
	const EdgeState* NetworkReliabilityObs::getState() const
	{
		return state.get();
	}
	NetworkReliabilityObs& NetworkReliabilityObs::operator=(const NetworkReliabilityObs& other)
	{
		if(&context != &other.context)
		{
			throw std::runtime_error("Internal error");
		}
		state = other.state;
		return *this;
	}
	NetworkReliabilityObs::NetworkReliabilityObs(NetworkReliabilityObs&& other)
		:context(other.context), state(other.state)
	{}
	const Context& NetworkReliabilityObs::getContext() const
	{
		return context;
	}
	void NetworkReliabilityObs::getReducedGraph(Context::internalGraph& outputGraph, int& minimumInoperative, std::vector<int>& edgeCounts, std::vector<int>& components, boost::detail::depth_first_visit_restricted_impl_helper<Context::internalGraph>::stackType& stack, std::vector<boost::default_color_type>& colorMap) const
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
	}
	void NetworkReliabilityObs::getReducedGraphNoSelfWithWeights(getReducedGraphNoSelfWithWeightsInput& input) const
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
			for(std::size_t i = 0; i < input.interestVertices.size(); i++)
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
	NetworkReliabilityObsWithContext::NetworkReliabilityObsWithContext(NetworkReliabilityObs& inputObs)
	{
		boost::shared_array<EdgeState> copiedState(new EdgeState[inputObs.getContext().getNEdges()]);
		memcpy(copiedState.get(), inputObs.getState(), sizeof(EdgeState) * inputObs.getContext().getNEdges());
		obs.reset(new NetworkReliabilityObs(inputObs.getContext(), copiedState));
	}
	const NetworkReliabilityObs& NetworkReliabilityObsWithContext::getObs() const
	{
		return *obs;
	}
	const Context& NetworkReliabilityObsWithContext::getContext() const
	{
		if(context) return *context;
		return obs->getContext();
	}

}
