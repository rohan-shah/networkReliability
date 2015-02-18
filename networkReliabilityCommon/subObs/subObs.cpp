#include "subObs/subObs.h"
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/math/distributions/binomial.hpp>
#include "graphAlgorithms.h"
#include "seriesParallelReduction.hpp"
namespace networkReliability
{
	namespace subObs
	{
		subObs::subObs(Context const& context, boost::shared_array<EdgeState> state, double radius)
			: ::networkReliability::NetworkReliabilityObs(context, state), radius(radius)
		{
		}
		void subObs::getReducedGraph(Context::internalGraph& outputGraph, std::vector<int>& edgeCounts, std::vector<int>& components, boost::detail::depth_first_visit_restricted_impl_helper<Context::internalGraph>::stackType& stack, std::vector<boost::default_color_type>& colorMap) const
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
		void subObs::getReducedGraphNoSelfWithWeights(getReducedGraphNoSelfWithWeightsInput& input) const
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
						input.reducedInterestVertices[i] = (int)*current;
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
		subObs::subObs(subObs&& other)
			: ::networkReliability::NetworkReliabilityObs(other.context, other.state), radius(other.radius)
		{
		}
		subObs& subObs::operator=(subObs&& other)
		{
			this->::networkReliability::NetworkReliabilityObs::operator=(static_cast<::networkReliability::NetworkReliabilityObs&&>(other));
			radius = other.radius;
			return *this;
		}
		double subObs::getRadius() const
		{
			return radius;
		}
	}
}
