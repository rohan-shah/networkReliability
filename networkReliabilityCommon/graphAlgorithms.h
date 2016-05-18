#ifndef GRAPH_ALGORITHMS_HEADER_GUARD
#define GRAPH_ALGORITHMS_HEADER_GUARD
#include "context.h"
#include "edgeState.h"
#include "depth_first_search_restricted.hpp"
namespace networkReliability
{
	typedef boost::color_traits<boost::default_color_type> Color;
	int countComponents(Context const& context, const edgeState* state, std::vector<int>& components, boost::detail::depth_first_visit_restricted_impl_helper<Context::internalGraph>::stackType& stack, std::vector<boost::default_color_type>& colorMap);
	bool isSingleComponent(Context const& context, const edgeState* state, std::vector<int>& components, boost::detail::depth_first_visit_restricted_impl_helper<Context::internalGraph>::stackType& stack, std::vector<boost::default_color_type>& colorMap);
	void identifyMinCutEdges(const std::vector<int>& edgeResidualCapacityVector, const std::vector<int>& capacityVector, const edgeState* state, const Context::internalDirectedGraph& directedGraph, std::vector<boost::default_color_type>& colorVector, std::vector<bool>& edgesToConsider, std::vector<int>& outputEdges, int source, int sink);
}
#endif
