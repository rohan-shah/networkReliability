#ifndef EXHAUSTIVE_SEARCH_HEADER_GUARD
#define EXHAUSTIVE_SEARCH_HEADER_GUARD
#include "context.h"
namespace networkReliability
{
	struct exhaustiveSearchArgs
	{
	public:
#ifdef _MSC_VER
		typedef long long counterType;
#else
		typedef unsigned long long counterType;
#endif
		exhaustiveSearchArgs(const context::internalGraph& graph)
			:graph(graph)
		{}
		const context::internalGraph graph;
		std::vector<counterType> result;
		std::vector<int> interestVertices;
		bool countDisconnected;
	};
	void exhaustiveSearch(exhaustiveSearchArgs& args);
}
#endif
