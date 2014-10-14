#include "graphAlgorithms.h"
#include "connected_components_restricted.hpp"
namespace networkReliability
{
	int countComponents(Context const& context, const EdgeState* state, std::vector<int>& components, boost::detail::depth_first_visit_restricted_impl_helper<Context::internalGraph>::stackType& stack, std::vector<boost::default_color_type>& colorMap)
	{
		const std::size_t nVertices = boost::num_vertices(context.getGraph());
		components.resize(nVertices);
		
		colorMap.clear();
		colorMap.resize(nVertices, Color::white());

		return boost::connected_components_restricted(context.getGraph(), &(components[0]), &(colorMap[0]), stack, state);
	}
	bool isSingleComponent(Context const& context, const EdgeState* state, std::vector<int>& components, boost::detail::depth_first_visit_restricted_impl_helper<Context::internalGraph>::stackType& stack, std::vector<boost::default_color_type>& colorMap)
	{
		countComponents(context, state, components, stack, colorMap);
		const std::vector<int>& interestVertices = context.getInterestVertices();
		int firstInterestComponentID = components[interestVertices[0]];
		for(int i = 1; i < interestVertices.size(); i++)
		{
			if(components[interestVertices[i]] != firstInterestComponentID)
			{
				return false;
			}
		}
		return true;
	}
}