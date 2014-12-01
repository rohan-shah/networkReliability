#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "Arguments.h"
#include "ArgumentsMPFR.h"
#include "graphAlgorithms.h"
#include <boost/graph/adjacency_matrix.hpp>
#include <iomanip>
#include <omp.h>
namespace networkReliability
{
	int main(int argc, char **argv)
	{
		omp_set_num_threads(6);

		boost::program_options::options_description options("Usage");
		options.add_options()
			("gridGraph", boost::program_options::value<int>(), "(int) The dimension of the square grid graph to use. Incompatible with graphFile and torusGraph. ")
			("graphFile", boost::program_options::value<std::string>(), "(string) The path to a graphml file. Incompatible with gridGraph")
			("completeGraph", boost::program_options::value<int>(), "(int) The number of vertices of the complete graph to use. ")
			("interestVertices", boost::program_options::value<std::vector<int> >()->multitoken(), "(int) The vertices of interest, that should be connected. ")
			("countDisconnected", boost::program_options::value<bool>()->default_value(false)->implicit_value(true), "(flag) Should we count the number of disconnected subgraphs?")
			("help", "Display this message");
		boost::program_options::variables_map variableMap;
		try
		{
			boost::program_options::store(boost::program_options::parse_command_line(argc, argv, options), variableMap);
		}
		catch(boost::program_options::error& ee)
		{
			std::cerr << "Error parsing command line arguments: " << ee.what() << std::endl << std::endl;
			std::cerr << options << std::endl;
			return -1;
		}

		if(variableMap.count("help") > 0)
		{
			std::cout << 
				"This program estimates the probability that the given graph is unreliable for the given vertices. That is, if edges are retained with a certain probability, what is the probability that the specified vertices are not all in the same connected component?\n\n"
			;
			std::cout << options << std::endl;
			return 0;
		}

		bool countDisconnected = variableMap["countDisconnected"].as<bool>();
		Context context = Context::emptyContext();
		if(!readContext(variableMap, context, 0.5))
		{
			return 0;
		}
		const std::size_t nEdges = context.getNEdges();

		if(nEdges > 36)
		{
			std::cout << "This program can be run with at most 36 edges" << std::endl;
			return 0;
		}
		Context::internalGraph const& graph = context.getGraph();
		const std::size_t nVertices = boost::num_vertices(graph);

		const std::vector<int> interestVertices = context.getInterestVertices();

		typedef long long counterType;
		const counterType maximumState = 1LL << nEdges;

		counterType* sizeCounters = new counterType[nEdges+1];
		memset(sizeCounters, 0, sizeof(counterType) * (nEdges+1));

		#pragma omp parallel
		{
			counterType* privateSizeCounters = new counterType[nEdges+1];
			memset(privateSizeCounters, 0, sizeof(counterType) * (nEdges+1));

			std::vector<int> privateConnectedComponents(nVertices);
			boost::detail::depth_first_visit_restricted_impl_helper<Context::internalGraph>::stackType stack;
			std::vector<boost::default_color_type> colorMap;

			std::vector<EdgeState> edgeStates(nEdges);
			EdgeState* edgeStatePtr = &(edgeStates[0]);

			#pragma omp for
			for(counterType state = 0; state < maximumState; state++)
			{
				int nEdgesThisGraph = 0;
				for(int edgeCounter = 0; edgeCounter < nEdges; edgeCounter++)
				{
					if(state & (1LL << edgeCounter))
					{
						edgeStatePtr[edgeCounter] = UNFIXED_OP;
						nEdgesThisGraph++;
					}
					else
					{
						edgeStatePtr[edgeCounter] = UNFIXED_INOP;
					}
				}
				bool currentGraphConnected = isSingleComponent(context, edgeStatePtr, privateConnectedComponents, stack, colorMap);
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
				for(int i = 0; i < nEdges+1; i++)
				{
					sizeCounters[i] += privateSizeCounters[i];
				}
			}
			delete[] privateSizeCounters;
		}
		
		//
		std::cout << "Command " << argv[0] << " run from directory \"" << boost::filesystem::current_path().string() << "\" with arguments \"";
		for(int i = 1; i < argc-1; i++)
		{
			std::cout << argv[i] << " "; 
		}
		std::cout << argv[argc-1] << "\"" << std::endl;

		if (countDisconnected) std::cout << "Number of disconnected subgraphs with that number of edges" << std::endl;
		else std::cout << "Number of connected subgraphs with that number of edges" << std::endl;
		for(int i = 0; i < nEdges+1; i++)
		{
			std::cout << std::setw(3) << i << ":  " << sizeCounters[i] << std::endl;
		}
		delete[] sizeCounters;
		return 0;
	}
}
int main(int argc, char **argv)
{
	return networkReliability::main(argc, argv);
}
