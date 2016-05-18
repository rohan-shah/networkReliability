#include <boost/program_options.hpp>
#include "arguments.h"
#include "argumentsMPFR.h"
#include <fstream>
#include "formulaDriver.h"
#include "empiricalDistribution.h"
#include <boost/archive/binary_iarchive.hpp>
namespace networkReliability
{
	int main(int argc, char **argv)
	{
		boost::program_options::options_description options("Usage");
		options.add_options()
			("distributionFile", boost::program_options::value<std::string>(), "(path) The path to an empirical distributions file. ")
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

		if(variableMap.count("distributionFile") != 1)
		{
			std::cout << "Please enter a single value for input `distributionFile'" << std::endl;
			return 0;
		}
		try
		{
			std::string distributionFile = variableMap["distributionFile"].as<std::string>();
			std::ifstream stream(distributionFile.c_str(), std::ios_base::binary);
			boost::archive::binary_iarchive archive(stream);
			empiricalDistribution loadedDistribution(archive);

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
			mpfr_class sumAllWeights = 0;
			for(std::size_t sampleCounter = 0; sampleCounter < sampleSize; sampleCounter++)
			{
				loadedDistribution.expand(sampleCounter, edges);
				double weight = 1;
				if(loadedDistribution.isWeighted()) weight = loadedDistribution.getWeight(sampleCounter);
				for(std::size_t i = 0; i < nFunctions; i++)
				{
					sums[i] += weight * driver.result[i]->calculate(edges);
				}
				sumAllWeights += weight;
			}
			for(std::size_t i = 0; i < nFunctions; i++)
			{
				std::cout << "Estimated value of function " << i << " was " << (sums[i]/sumAllWeights) << std::endl;
			}
		}
		catch(std::runtime_error& err)
		{
			std::cout << "Runtime error: " << err.what() << std::endl;
			return 0;
		}
		return 0;
	}
}
int main(int argc, char **argv)
{
	return networkReliability::main(argc, argv);
}
