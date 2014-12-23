#include <boost/program_options.hpp>
#include "Arguments.h"
#include "ArgumentsMPFR.h"
#include <fstream>
#include "formulaDriver.h"
#include "empiricalDistribution.h"
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
		empiricalDistribution loadedDistribution = empiricalDistribution::load(distributionFile);
		const std::size_t nEdges = loadedDistribution.getNEdges();
		const std::size_t sampleSize = loadedDistribution.getNSamples();

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
		
		for(std::size_t sampleCounter = 0; sampleCounter < sampleSize; sampleCounter++)
		{
			loadedDistribution.expand(sampleCounter, edges);
			for(int i = 0; i < nFunctions; i++)
			{
				sums[i] += driver.result[i]->calculate(edges);
			}
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
