#include "crudeMCImpl.h"
#include "networkReliabilityObs.h"
#include "graphAlgorithms.h"
namespace networkReliability
{
	double crudeMC(crudeMCArgs& args)
	{
		std::size_t n = args.n;
		if(n <= 0)
		{
			throw std::runtime_error("Input n must be positive");
		}
		
		const context& contextObj = args.contextObj;
		boost::detail::depth_first_visit_restricted_impl_helper<context::internalGraph>::stackType stack;

		std::vector<int> components;
		std::vector<int> interestComponents;
		const std::vector<int> interestVertices = contextObj.getInterestVertices();

		std::vector<boost::default_color_type> colorMap;
		int countDisconnected = 0;

		for(std::size_t i = 0; i < n; i++)
		{
			NetworkReliabilityObs obs(contextObj, args.randomSource);
			if(!isSingleComponent(contextObj.getGraph(), obs.getState(), components, stack, colorMap, contextObj.getInterestVertices()))
			{
				countDisconnected++;
			}			
		}
		return (float)countDisconnected / (float)n;
	}
}
