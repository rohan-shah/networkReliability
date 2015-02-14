#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "Arguments.h"
#include "ArgumentsMPFR.h"
#include "graphAlgorithms.h"
#include <boost/graph/adjacency_matrix.hpp>
#include <iomanip>
#include <omp.h>
#include "formulaDriver.h"
#include <fstream>
#include <dlfcn.h>
#include "createFunctionBinary.h"
#include <boost/ptr_container/ptr_vector.hpp>
#include "Context.h"
namespace networkReliability
{
	int main(int argc, char **argv)
	{
		boost::program_options::options_description options("Usage");
		options.add_options()
			("gridGraph", boost::program_options::value<int>(), "(int) The dimension of the square grid graph to use. Incompatible with graphFile and torusGraph. ")
			("graphFile", boost::program_options::value<std::string>(), "(string) The path to a graphml file. Incompatible with gridGraph")
			("completeGraph", boost::program_options::value<int>(), "(int) The number of vertices of the complete graph to use. ")
			("interestVertices", boost::program_options::value<std::vector<int> >()->multitoken(), "(int) The vertices of interest, that should be connected. ")
			("function", boost::program_options::value<std::string>(), "(string) The function to evaluate.")
			("functionFile", boost::program_options::value<std::string>(), "(string) File containing the function to evaluate")
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
		
		std::string message, function;
		formulaDriver driver(nEdges);
		if(!readFunctionFile(variableMap, driver, function, message))
		{
			std::cout << message << std::endl;
			return 0;
		}
		std::size_t nFunctions = driver.result.size();
		//Work out the range of the function
		std::vector<int> functionMax(nFunctions, 0), functionMin(nFunctions, 0), functionRanges(nFunctions, 0);
		for(std::size_t i = 0; i < nFunctions; i++)
		{
			int& currentFunctionMax = functionMax[i];
			int& currentFunctionMin = functionMin[i];
			std::vector<int>& currentEdgeIDs = driver.edgeIDs[i];
			std::vector<int> edgeValues(nEdges, 0);
			boost::shared_ptr<node> currentResult = driver.result[i];
			currentFunctionMax = currentFunctionMin = currentResult->calculate(edgeValues);
			{
				long maxMaskValue = (1 << currentEdgeIDs.size());
				for(long maskCounter = 0; maskCounter < maxMaskValue; maskCounter++)
				{
					for(std::vector<int>::iterator relevantEdgesIterator = currentEdgeIDs.begin(); relevantEdgesIterator != currentEdgeIDs.end(); relevantEdgesIterator++)
					{
						int currentBit = std::distance(currentEdgeIDs.begin(), relevantEdgesIterator);
						edgeValues[*relevantEdgesIterator] = (maskCounter & (1 << currentBit)) >> currentBit;
					}
					int value = currentResult->calculate(edgeValues);
					currentFunctionMin = std::min(currentFunctionMin, value);
					currentFunctionMax = std::max(currentFunctionMax, value);
				}
			}
			functionRanges[i] = currentFunctionMax - currentFunctionMin+1;
		}
		std::string binaryFile;
		if(!createFunctionsBinary(function, driver.edgeIDs, message, binaryFile))
		{
			std::cout << "Error calling createFunctionBinary: " << message << std::endl;
			return 0;
		}
		void* compiledHandle = dlopen(binaryFile.c_str(), RTLD_LAZY);
		if(compiledHandle == NULL)
		{
			std::cout << "Unable to load compiled function code:" << dlerror() << std::endl;
			return 0;
		}
		typedef int (*compiledFunctionType)(long);
		std::vector<compiledFunctionType> compiledFunctions(nFunctions, NULL);
		for(std::size_t i = 0; i < nFunctions; i++)
		{
			std::string functionName = "compiledFunction" + boost::lexical_cast<std::string>(i);
			compiledFunctionType compiledFunction = (compiledFunctionType)dlsym(compiledHandle, functionName.c_str());
			if(compiledFunction == NULL)
			{
				std::cout << "Unable to resolve symbol " << functionName << " in generated code" << std::endl;
				return 0;
			}
			compiledFunctions[i] = compiledFunction;
		}
		//safe to unlink already, it'll stay in memory.
		unlink(binaryFile.c_str());

		typedef long long counterType;
		const counterType maximumState = 1LL << nEdges;

		boost::ptr_vector<counterType> sizeCounterArrays;
		for(std::vector<int>::iterator functionRangesIterator = functionRanges.begin(); functionRangesIterator != functionRanges.end(); functionRangesIterator++) 
		{
			counterType* currentArray = new counterType[(nEdges+1) * *functionRangesIterator];
			sizeCounterArrays.push_back(currentArray);
			memset(currentArray, 0, sizeof(counterType) * (nEdges+1)* *functionRangesIterator);
		}
		#pragma omp parallel
		{
			boost::ptr_vector<counterType> privateSizeCounterArrays;
			for(std::vector<int>::iterator functionRangesIterator = functionRanges.begin(); functionRangesIterator != functionRanges.end(); functionRangesIterator++)
			{
				counterType* currentArray = new counterType[(nEdges+1)* *functionRangesIterator];
				memset(currentArray, 0, sizeof(counterType) * (nEdges+1)* *functionRangesIterator);
				privateSizeCounterArrays.push_back(currentArray);
			}

			std::vector<int> privateConnectedComponents(nVertices);
			boost::detail::depth_first_visit_restricted_impl_helper<Context::internalGraph>::stackType stack;
			std::vector<boost::default_color_type> colorMap;

			std::vector<EdgeState> edgeStates(nEdges);
			EdgeState* edgeStatePtr = &(edgeStates[0]);

			#pragma omp for
			for(counterType state = 0; state < maximumState; state++)
			{
				int nEdgesThisGraph = 0;
				for(std::size_t edgeCounter = 0; edgeCounter < nEdges; edgeCounter++)
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
				if(!currentGraphConnected)
				{
					for(std::size_t i = 0; i < nFunctions; i++)
					{
						int functionValue = compiledFunctions[i](state);
						(&(privateSizeCounterArrays[i]))[nEdgesThisGraph*functionRanges[i] + (functionValue - functionMin[i])]++;
					}
				}
			}
			#pragma omp critical
			{
				for(std::size_t j = 0; j < nFunctions; j++)
				{
					for(std::size_t i = 0; i < (nEdges+1)*functionRanges[j]; i++)
					{
						(&(sizeCounterArrays[j]))[i] += (&(privateSizeCounterArrays[j]))[i];
					}
				}
			}
		}
		
		//
		std::cout << "Command " << argv[0] << " run from directory \"" << boost::filesystem::current_path().string() << "\" with arguments \"";
		for(int i = 1; i < argc-1; i++)
		{
			std::cout << argv[i] << " "; 
		}
		std::cout << argv[argc-1] << "\"" << std::endl;
		std::cout << nFunctions << " functions input" << std::endl;
		std::cout << "Graph had " << nEdges << " edges" << std::endl;
		for(std::size_t j = 0; j < nFunctions;j ++)
		{
			std::cout << "Function " << j << " took on " << functionRanges[j] << " values" << std::endl;
			for(std::size_t i = 0; i < (nEdges+1)*functionRanges[j]; i++)
			{
				int functionValue = (i % functionRanges[j]);
				int nEdges = i / functionRanges[j];
				std::cout << std::setw(3) << nEdges << " edges, value " << std::setw(3) << (functionValue+functionMin[j])  << ":  " << (&(sizeCounterArrays[j]))[nEdges*functionRanges[j] + functionValue] << std::endl;
			}
		}
		return 0;
	}
}
int main(int argc, char **argv)
{
	return networkReliability::main(argc, argv);
}
