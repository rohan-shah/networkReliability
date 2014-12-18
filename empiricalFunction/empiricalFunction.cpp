#include <boost/program_options.hpp>
#include "Arguments.h"
#include "ArgumentsMPFR.h"
#include <fstream>
#include "formulaDriver.h"
namespace networkReliability
{
	int main(int argc, char **argv)
	{
		boost::program_options::options_description options("Usage");
		options.add_options()
			("distributionFile", boost::program_options::value<std::string>(), "(path) The path to an empirical distributions file. ")
			("opProbability", boost::program_options::value<std::string>(), "(float) The probability that an edge is operational. ")
			("function", boost::program_options::value<std::string>(), "(string) The function to evaluate.")
			("functionFile", boost::program_options::value<std::string>(), "(string) File containing the function to evaluate");

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

		mpfr_class opProbability;
		if(!readProbabilityString(variableMap, opProbability))
		{
			std::cout << "Unable to read input `opProbability'" << std::endl;
			return 0;
		}
		if(variableMap.count("distributionFile") != 1)
		{
			std::cout << "Please enter a single value for input `distributionFile'" << std::endl;
			return 0;
		}
		std::string distributionFile = variableMap["distributionFile"].as<std::string>();
		std::ifstream distributionsStream(distributionFile, std::ios_base::in | std::ios_base::binary);
		if(!distributionsStream)
		{
			std::cout << "Unable to open specified file" << std::endl;
			return 0;
		}
		std::size_t nEdges, sampleSize;
		distributionsStream.read((char*)&nEdges, sizeof(std::size_t));
		distributionsStream.read((char*)&sampleSize, sizeof(std::size_t));

		std::string function, message;
		formulaDriver driver(nEdges);
		
		if(!readFunctionFile(variableMap, driver, function, message))
		{
			std::cout << message << std::endl;
			return 0;
		}
		const std::size_t nFunctions = driver.result.size();
		//States of individual edges
		std::vector<int> edges(nEdges, 0);
		std::vector<mpfr_class> sums(nEdges, 0);
		
		//read by bytes
		/*for(int j = 0;j < sampleSize; j++)
		{
			for(int i = 0; i < nEdges; i++)
			{
				char c;
				distributionsStream >> c;
				edges[i] = c;
			}
			for(int i = 0; i < nFunctions; i++)
			{
				sums[i] += driver.result[i]->calculate(edges);
			}
		}*/
		//read by bits
		int nStoredBits = 0;
		unsigned int currentValue = 0;
		for(int j = 0;j < sampleSize; j++)
		{
			for(int i = 0; i < nEdges; i++)
			{
				if(!nStoredBits)
				{
					distributionsStream.read((char*)&currentValue, sizeof(int));
					nStoredBits = sizeof(int)*8;
				}
				edges[i] = (currentValue & (1U << (8*sizeof(int)-1))) > 0;
				currentValue<<=1;
				nStoredBits--;
			}
			for(int i = 0; i < nFunctions; i++)
			{
				sums[i] += driver.result[i]->calculate(edges);
			}
		}
		int position = distributionsStream.tellg();
		distributionsStream.seekg(0, distributionsStream.end);
		int length = distributionsStream.tellg();
		if(position != length)
		{
			std::cout << "Empirical distributions file was too long. Finished reading at position " << position << " of " << length << std::endl;
			return 0;
		}
		for(int i = 0; i < nFunctions; i++)
		{
			std::cout << "Estimated value of function " << i << " was " << (sums[i]/sampleSize) << std::endl;
		}
		return 0;
	}
}
int main(int argc, char **argv)
{
	return networkReliability::main(argc, argv);
}
