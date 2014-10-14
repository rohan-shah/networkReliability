#ifndef GRAPH_ALGORITHMS_HEADER_GUARD
#define GRAPH_ALGORITHMS_HEADER_GUARD
#include "Context.h"
#include "EdgeState.h"
#include "depth_first_search_restricted.hpp"
namespace networkReliability
{
	typedef boost::color_traits<boost::default_color_type> Color;
	int countComponents(Context const& context, const EdgeState* state, std::vector<int>& components, boost::detail::depth_first_visit_restricted_impl_helper<Context::internalGraph>::stackType& stack, std::vector<boost::default_color_type>& colorMap);
	bool isSingleComponent(Context const& context, const EdgeState* state, std::vector<int>& components, boost::detail::depth_first_visit_restricted_impl_helper<Context::internalGraph>::stackType& stack, std::vector<boost::default_color_type>& colorMap);
}
#endif