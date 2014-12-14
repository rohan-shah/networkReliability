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
namespace networkReliability
{
	std::string createFunctionBinary(std::string function, const std::vector<int>& edgeIDs)
	{
		char codeFile[] = "./codeFileXXXXXX";
		mkstemps(codeFile, 0);
		std::ofstream outputCodeStream(codeFile, std::ios::out);
		outputCodeStream << "int compiledFunction(long edgeMask)" << std::endl;
		outputCodeStream << "{" << std::endl;
		for(std::vector<int>::const_iterator relevantEdgeIterator = edgeIDs.begin(); relevantEdgeIterator != edgeIDs.end(); relevantEdgeIterator++)
		{
			outputCodeStream << "\tint e" << *relevantEdgeIterator << " = (edgeMask & (1 << " << *relevantEdgeIterator << ")) >> " << *relevantEdgeIterator << ";" << std::endl;
		}
		outputCodeStream << "\treturn " << function << ";" << std::endl;
		outputCodeStream << "}" << std::endl;
		outputCodeStream.flush();
		outputCodeStream.close();
		char binaryFile[] = "./binaryFileXXXXXX.o";
		mkstemps(binaryFile, 2);

		std::string binaryFileString = binaryFile;
		//Compile .o file
		std::stringstream commandStream;
		commandStream << "gcc -fpic -x c -O3 -c " << codeFile << " -o " << binaryFileString;
		system(commandStream.str().c_str());
		commandStream.str("");
		commandStream.clear();
		//We don't need code file any more
		unlink(codeFile);

		char sharedObjectFile[] = "./binaryFileXXXXXX.so";
		mkstemps(sharedObjectFile, 3);
		std::string sharedObjectString = sharedObjectFile;
		//compile .so file
		commandStream << "gcc -shared -o " << sharedObjectString << " "  << binaryFileString;
		system(commandStream.str().c_str());
		//We don't need .o file any more
		unlink(binaryFileString.c_str());
		return sharedObjectString;
	}
	int main(int argc, char **argv)
	{
		omp_set_num_threads(6);

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
		if(variableMap.count("function") + variableMap.count("functionFile") != 1)
		{
			std::cout << "Exactly one of `function' or `functionFile' is requied" << std::endl;
		}

		if(nEdges > 36)
		{
			std::cout << "This program can be run with at most 36 edges" << std::endl;
			return 0;
		}
		Context::internalGraph const& graph = context.getGraph();
		const std::size_t nVertices = boost::num_vertices(graph);
		const std::vector<int> interestVertices = context.getInterestVertices();
		
		std::string functionFile, function;
		if(variableMap.count("function") > 0)
		{
			function = variableMap["function"].as<std::string>();
			char tmpFunctionFile[] = "./functionFileXXXXXX";
			mkstemp(tmpFunctionFile);
			functionFile = tmpFunctionFile;
			FILE* tmpFunctionFileHandle = fopen(functionFile.c_str(), "w");
			fwrite(function.c_str(), sizeof(char), function.size(), tmpFunctionFileHandle);
			fclose(tmpFunctionFileHandle);
		}
		else 
		{
			functionFile = variableMap["functionFile"].as<std::string>();
			std::ifstream functionFileStream(functionFile, std::ios::in);
			if(!functionFileStream)
			{
				std::cout << "Unable to read data from function file" << std::endl;
				return 0;
			}
			function = std::string(std::istreambuf_iterator<char>(functionFileStream), std::istreambuf_iterator<char>());
		}

		formulaDriver driver(nEdges);
		int parseResult = driver.parse(functionFile);
		//Unlike the temporary file, if we used one. 
		if(variableMap.count("function")) unlink(functionFile.c_str());
		//Get out the unique edges that are of interest for this function. 
		std::sort(driver.edgeIDs.begin(), driver.edgeIDs.end());
		driver.edgeIDs.erase(std::unique(driver.edgeIDs.begin(), driver.edgeIDs.end()), driver.edgeIDs.end());
		if(parseResult != 0 || driver.message.size() != 0)
		{
			std::cout << "Errors parsing function:" << std::endl;
			std::cout << driver.message;
			return 0;
		}
		//Work out the range of the function
		std::vector<int> edgeValues(nEdges, 0);
		int functionMax, functionMin;
		functionMax = functionMin = driver.result->calculate(edgeValues);
		{
			long maxMaskValue = (1 << driver.edgeIDs.size());
			for(long maskCounter = 0; maskCounter < maxMaskValue; maskCounter++)
			{
				for(std::vector<int>::iterator relevantEdgesIterator = driver.edgeIDs.begin(); relevantEdgesIterator != driver.edgeIDs.end(); relevantEdgesIterator++)
				{
					int currentBit = std::distance(driver.edgeIDs.begin(), relevantEdgesIterator);
					edgeValues[*relevantEdgesIterator] = (maskCounter & (1 << currentBit)) >> currentBit;
				}
				int value = driver.result->calculate(edgeValues);
				functionMin = std::min(functionMin, value);
				functionMax = std::max(functionMax, value);
			}
		}
		const int functionValues = functionMax - functionMin + 1;
		//If it's constant valued, return immediately
		if(functionMin == functionMax)
		{
			std::cout << "Expected value was " << functionMax << std::endl;
			return 0;
		}
		std::string binaryFile = createFunctionBinary(function, driver.edgeIDs);
		void* compiledHandle = dlopen(binaryFile.c_str(), RTLD_LAZY);
		if(compiledHandle == NULL)
		{
			std::cout << "Unable to load compiled function code:" << dlerror() << std::endl;
			return 0;
		}
		typedef int (*compiledFunctionType)(long);
		compiledFunctionType compiledFunction = (compiledFunctionType)dlsym(compiledHandle, "compiledFunction");
		if(compiledFunction == NULL)
		{
			std::cout << "Unable to resolve symbol compiledFunction in generated code" << std::endl;
			return 0;
		}
		//safe to unlink already, it'll stay in memory.
		unlink(binaryFile.c_str());

		typedef long long counterType;
		const counterType maximumState = 1LL << nEdges;

		counterType* sizeCounters = new counterType[(nEdges+1) * functionValues];
		memset(sizeCounters, 0, sizeof(counterType) * (nEdges+1)*functionValues);

		#pragma omp parallel
		{
			counterType* privateSizeCounters = new counterType[(nEdges+1)*functionValues];
			memset(privateSizeCounters, 0, sizeof(counterType) * (nEdges+1)*functionValues);

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
				if(currentGraphConnected)
				{
					int functionValue = compiledFunction(state);
					privateSizeCounters[nEdgesThisGraph*functionValues + (functionValue - functionMin)]++;
				}
			}
			#pragma omp critical
			{
				for(int i = 0; i < (nEdges+1)*functionValues; i++)
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
		std::cout << "Graph had " << nEdges << " edges" << std::endl;
		std::cout << "Function took on " << functionValues << " values" << std::endl;

		for(int i = 0; i < (nEdges+1)*functionValues; i++)
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
