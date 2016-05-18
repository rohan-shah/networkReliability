#include "graphAlgorithms.h"
#include "connected_components_restricted.hpp"
#include <boost/graph/filtered_graph.hpp>
namespace networkReliability
{
	int countComponents(Context const& context, const edgeState* state, std::vector<int>& components, boost::detail::depth_first_visit_restricted_impl_helper<Context::internalGraph>::stackType& stack, std::vector<boost::default_color_type>& colorMap)
	{
		const std::size_t nVertices = boost::num_vertices(context.getGraph());
		components.resize(nVertices);
		
		colorMap.clear();
		colorMap.resize(nVertices, Color::white());

		return boost::connected_components_restricted(context.getGraph(), &(components[0]), &(colorMap[0]), stack, state);
	}
	bool isSingleComponent(Context const& context, const edgeState* state, std::vector<int>& components, boost::detail::depth_first_visit_restricted_impl_helper<Context::internalGraph>::stackType& stack, std::vector<boost::default_color_type>& colorMap)
	{
		countComponents(context, state, components, stack, colorMap);
		const std::vector<int>& interestVertices = context.getInterestVertices();
		int firstInterestComponentID = components[interestVertices[0]];
		for(std::size_t i = 1; i < interestVertices.size(); i++)
		{
			if(components[interestVertices[i]] != firstInterestComponentID)
			{
				return false;
			}
		}
		return true;
	}
	class identify_dfs_visitor : public boost::default_dfs_visitor
	{
	public:
		struct linkedListTerm
		{
			int parentIndex;
			int value;
		};
		identify_dfs_visitor(const std::vector<int>& edgeResidualCapacityVector, const std::vector<int>& edgeCapacityVector, const Context::internalDirectedGraph& directedGraph, const edgeState* state, std::vector<linkedListTerm>& perVertex)
			: edgeResidualCapacityVector(edgeResidualCapacityVector), state(state), edgeCapacityVector(edgeCapacityVector), perVertex(perVertex)
		{
			std::size_t nVertices = boost::num_vertices(directedGraph);
			perVertex.resize(nVertices);
			for(std::size_t i = 0; i < nVertices; i++)
			{
				perVertex[i].value = INT_MIN;
				perVertex[i].parentIndex = -1;
			}
		}
		void join(linkedListTerm* first, linkedListTerm* second)
		{
			linkedListTerm newTerm;
			newTerm.parentIndex = -1;
			newTerm.value = second->value;
			if(newTerm.value == INT_MIN) newTerm.value = first->value;
			first->parentIndex = second->parentIndex = perVertex.size();
			first->value = second->value = -1;
			perVertex.push_back(newTerm);
		}
		template<typename Edge, typename Graph> void tree_edge(Edge e, Graph& g)
		{
			doEdge(e, g);
		}
		template<typename Edge, typename Graph> void doEdge(Edge e, Graph& g)
		{
			int index = boost::get(boost::edge_index, g, e);
			int reverseIndex;
			if(index % 2) reverseIndex = index-1;
			else reverseIndex = index+1;
			int edgeResidualFlow = edgeResidualCapacityVector[index];
			if(edgeResidualFlow > edgeCapacityVector[index]) edgeResidualFlow = edgeCapacityVector[index];
			int edgeReverseResidualFlow = edgeResidualCapacityVector[reverseIndex];
			if(edgeReverseResidualFlow > edgeCapacityVector[reverseIndex]) edgeReverseResidualFlow = edgeCapacityVector[reverseIndex];

			int edgeFlow = edgeCapacityVector[index] - edgeResidualFlow;
			int edgeReverseFlow = edgeCapacityVector[reverseIndex] - edgeReverseResidualFlow;
			linkedListTerm* first = getTerminal(boost::source(e, g));
			linkedListTerm* second = getTerminal(boost::target(e, g));
			//If there's no flow either direction, both vertices should have the same label
			if((edgeFlow == 0 && edgeReverseFlow == 0) || state[index/2] == FIXED_OP)
			{
				if(first != second)
				{
					if(first->value == INT_MIN && second->value == INT_MIN) first->value = 0;
					join(first, second);
				}
			}
			//If we're going AGAINST the flow, it should decrease by one (at most)
			else if(edgeFlow > 0 && edgeReverseFlow == 0)
			{
				if(first != second)
				{
					if(first->value == INT_MIN && second->value == INT_MIN)
					{
						first->value = 0;
						second->value = -1;
					}
					if(first->value == INT_MIN) first->value = second->value+1;
					if(second->value == INT_MIN) second->value = first->value -1;
					if(second->value < first->value - 1)
					{
						linkedListTerm newTerm;
						newTerm.parentIndex = -1;
						newTerm.value = first->value - 1;
						second->parentIndex = perVertex.size();
						second->value = -1;
						perVertex.push_back(newTerm);
					}
				}
			}
			//this corresponds to a case where flows swap over unneccessarily, so it counts as a zero flow edge
			else if(edgeFlow > 0 && edgeReverseFlow > 0)
			{
				if(first != second)
				{
					if(first->value == INT_MIN && second->value == INT_MIN)
					{
						first->value = 0;
						second->value = 0;
					}
					if(first->value == INT_MIN) first->value = second->value;
					if(second->value == INT_MIN) second->value = first->value;
					join(first, second);
				}
			}
			//This corresponds to going against the flow
			else
			{
				if(first->value >= second->value)
				{
					second->value = first->value;
					join(first, second);
				}
			}
		}
		template<typename Edge, typename Graph> void back_edge(Edge e, Graph& g)
		{
			doEdge(e, g);
		}
		template<typename Edge, typename Graph> void forward_or_cross_edge(Edge e, Graph& g)
		{
			doEdge(e, g);
		}
		linkedListTerm* getTerminal(int index)
		{
			linkedListTerm* ret = &(perVertex[index]);
			while(ret->parentIndex != -1) ret = &(perVertex[ret->parentIndex]);
			return ret;
		}
		const std::vector<int>& edgeResidualCapacityVector;
		const edgeState* state;
		const std::vector<int>& edgeCapacityVector;
		std::vector<linkedListTerm>& perVertex;
	};
	struct withFlowFilter
	{
	public:
		withFlowFilter(const Context::internalDirectedGraph& directedGraph, const std::vector<int>& edgeCapacityVector, const std::vector<int>& edgeResidualCapacityVector, const edgeState* state)
			:directedGraph(&directedGraph), edgeCapacityVector(&edgeCapacityVector), edgeResidualCapacityVector(&edgeResidualCapacityVector), state(state)
		{}
		template <typename Edge> bool operator()(const Edge& e) const
		{
			int index = boost::get(boost::edge_index, *directedGraph, e);
			if(state[index/2] == FIXED_INOP) return false;
			int reverseIndex;
			if(index % 2) reverseIndex = index-1;
			else reverseIndex = index+1;

			int edgeResidualFlow = (*edgeResidualCapacityVector)[index];
			if(edgeResidualFlow > (*edgeCapacityVector)[index]) edgeResidualFlow = (*edgeCapacityVector)[index];

			int edgeReverseResidualFlow = (*edgeResidualCapacityVector)[reverseIndex];
			if(edgeReverseResidualFlow > (*edgeCapacityVector)[reverseIndex]) edgeReverseResidualFlow = (*edgeCapacityVector)[reverseIndex];
			int edgeReverseFlow = (*edgeCapacityVector)[reverseIndex] - edgeReverseResidualFlow;
			int edgeFlow = (*edgeCapacityVector)[index] - edgeResidualFlow;
			//If flow goes in both directions and it's not one of the edges marked as HIGH_CAPACITY, then really it counts as an edge with zero flow in both directions
			if((edgeResidualFlow < HIGH_CAPACITY && edgeReverseResidualFlow < HIGH_CAPACITY) && edgeFlow > 0 && edgeReverseFlow > 0) return true;

			return edgeReverseFlow == 0 || edgeResidualFlow == HIGH_CAPACITY || edgeReverseResidualFlow == HIGH_CAPACITY;
		}
		withFlowFilter()
		{}
		const Context::internalDirectedGraph* directedGraph;
		const std::vector<int>* edgeCapacityVector;
		const std::vector<int>* edgeResidualCapacityVector;
		const edgeState* state;
	};
	struct againstFlowFilter
	{
	public:
		againstFlowFilter(const Context::internalDirectedGraph& directedGraph, const std::vector<int>& edgeCapacityVector, const std::vector<int>& edgeResidualCapacityVector, const edgeState* state)
			:directedGraph(&directedGraph), edgeCapacityVector(&edgeCapacityVector), edgeResidualCapacityVector(&edgeResidualCapacityVector), state(state)
		{}
		template <typename Edge> bool operator()(const Edge& e) const
		{
			int index = boost::get(boost::edge_index, *directedGraph, e);
			if(state[index/2] == FIXED_INOP) return false;
			int reverseIndex;
			if(index % 2) reverseIndex = index-1;
			else reverseIndex = index+1;

			int edgeResidualFlow = (*edgeResidualCapacityVector)[index];
			if(edgeResidualFlow > (*edgeCapacityVector)[index]) edgeResidualFlow = (*edgeCapacityVector)[index];

			int edgeReverseResidualFlow = (*edgeResidualCapacityVector)[reverseIndex];
			if(edgeReverseResidualFlow > (*edgeCapacityVector)[reverseIndex]) edgeReverseResidualFlow = (*edgeCapacityVector)[reverseIndex];
			int edgeReverseFlow = (*edgeCapacityVector)[reverseIndex] - edgeReverseResidualFlow;
			int edgeFlow = (*edgeCapacityVector)[index] - edgeResidualFlow;
			//If flow goes in both directions and it's not one of the edges marked as HIGH_CAPACITY, then really it counts as an edge with zero flow in both directions
			if((edgeResidualFlow < HIGH_CAPACITY && edgeReverseResidualFlow < HIGH_CAPACITY) && edgeFlow > 0 && edgeReverseFlow > 0) return true;
			return edgeFlow == 0 || edgeResidualFlow == HIGH_CAPACITY || edgeReverseResidualFlow == HIGH_CAPACITY;
		}
		againstFlowFilter()
		{}
		const Context::internalDirectedGraph* directedGraph;
		const std::vector<int>* edgeCapacityVector;
		const std::vector<int>* edgeResidualCapacityVector;
		const edgeState* state;

	};
	void identifyMinCutEdges(const std::vector<int>& edgeResidualCapacityVector, const std::vector<int>& capacityVector, const edgeState* state, const Context::internalDirectedGraph& directedGraph, std::vector<boost::default_color_type>& colorVector, std::vector<bool>& edgesToConsider, std::vector<int>& outputEdges, int source, int sink)
	{
		typedef boost::property_map<Context::internalDirectedGraph, boost::vertex_index_t>::const_type vertexIndexMapType;
		typedef boost::iterator_property_map<typename std::vector<boost::default_color_type>::iterator, vertexIndexMapType> colorMapType;
		vertexIndexMapType vertexIndexMap = boost::get(boost::vertex_index, directedGraph);
		colorMapType colorMap(colorVector.begin(), vertexIndexMap);

		std::vector<identify_dfs_visitor::linkedListTerm> perVertex;

		withFlowFilter withFlowFilterObj(directedGraph, capacityVector, edgeResidualCapacityVector, state);
		boost::filtered_graph<Context::internalDirectedGraph, withFlowFilter> withFlowFilteredGraph(directedGraph, withFlowFilterObj);
		identify_dfs_visitor vis(edgeResidualCapacityVector, capacityVector, directedGraph, state, perVertex);
		boost::depth_first_search(withFlowFilteredGraph, boost::visitor(vis).color_map(colorMap).root_vertex(source));

		againstFlowFilter againstFlowFilterObj(directedGraph, capacityVector, edgeResidualCapacityVector, state);
		boost::filtered_graph<Context::internalDirectedGraph, againstFlowFilter> againstFlowFilteredGraph(directedGraph, againstFlowFilterObj);
		boost::depth_first_search(againstFlowFilteredGraph, boost::visitor(vis).color_map(colorMap).root_vertex(sink));

		Context::internalDirectedGraph::edge_iterator current, end;
		boost::tie(current, end) = boost::edges(directedGraph);
		for(;current != end; current++)
		{
			int currentEdgeIndex = boost::get(boost::edge_index, directedGraph, *current);
			if(currentEdgeIndex % 2 && edgesToConsider[currentEdgeIndex/2])
			{
				identify_dfs_visitor::linkedListTerm* sourceTerm = vis.getTerminal(boost::source(*current, directedGraph));
				identify_dfs_visitor::linkedListTerm* targetTerm = vis.getTerminal(boost::target(*current, directedGraph));
				if(targetTerm->value != sourceTerm->value) outputEdges.push_back(currentEdgeIndex/2);
			}
		}
	}
}
