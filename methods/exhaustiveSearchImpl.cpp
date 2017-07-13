#include "exhaustiveSearchImpl.h"
#include "graphAlgorithms.h"
#include <boost/graph/adjacency_matrix.hpp>
#include <omp.h>
namespace networkReliability
{
	void exhaustiveSearch(exhaustiveSearchArgs& args)
	{
		bool countDisconnected = args.countDisconnected;
		const std::size_t nEdges = boost::num_edges(args.graph);

		if(nEdges > 36)
		{
			throw std::runtime_error("This program can be run with at most 36 edges");
		}
		context::internalGraph const& graph = args.graph;
		const std::size_t nVertices = boost::num_vertices(graph);

		const std::vector<int> interestVertices = args.interestVertices;
		for(std::vector<int>::const_iterator i = interestVertices.begin(); i != interestVertices.end(); i++)
		{
			if(*i < 0) throw std::runtime_error("Interest vertices must be non-negative");
			if(*i >= (int)nVertices) throw std::runtime_error("Input vertex was too large");
		}

		const exhaustiveSearchArgs::counterType maximumState = ((exhaustiveSearchArgs::counterType)1) << nEdges;

		args.result.resize(nEdges+1, 0);
		#pragma omp parallel
		{
			exhaustiveSearchArgs::counterType* privateSizeCounters = new exhaustiveSearchArgs::counterType[nEdges+1];
			memset(privateSizeCounters, 0, sizeof(exhaustiveSearchArgs::counterType) * (nEdges+1));

			std::vector<int> privateConnectedComponents(nVertices);
			boost::detail::depth_first_visit_restricted_impl_helper<context::internalGraph>::stackType stack;
			std::vector<boost::default_color_type> colorMap;

			std::vector<edgeState> edgeStates(nEdges);
			edgeState* edgeStatePtr = &(edgeStates[0]);

			#pragma omp for
			for(exhaustiveSearchArgs::counterType state = 0; state < maximumState; state++)
			{
				int nEdgesThisGraph = 0;
				for(std::size_t edgeCounter = 0; edgeCounter < nEdges; edgeCounter++)
				{
					if(state & (((exhaustiveSearchArgs::counterType)1) << edgeCounter))
					{
						edgeStatePtr[edgeCounter] = UNFIXED_OP;
						nEdgesThisGraph++;
					}
					else
					{
						edgeStatePtr[edgeCounter] = UNFIXED_INOP;
					}
				}
				bool currentGraphConnected = isSingleComponent(graph, edgeStatePtr, privateConnectedComponents, stack, colorMap, interestVertices) ;
				if(!countDisconnected && currentGraphConnected)
				{
					privateSizeCounters[nEdgesThisGraph]++;
				}
				else if (countDisconnected && !currentGraphConnected)
				{
					privateSizeCounters[nEdgesThisGraph]++;
				}
			}
			#pragma omp critical
			{
				for(std::size_t i = 0; i < nEdges+1; i++)
				{
					args.result[i] += privateSizeCounters[i];
				}
			}
			delete[] privateSizeCounters;
		}
	}
}
