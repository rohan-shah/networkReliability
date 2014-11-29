#include "Arguments.h"
#include <boost/iterator/counting_iterator.hpp>
namespace networkReliability
{
	bool readN(boost::program_options::variables_map& variableMap, std::size_t& out)
	{
		if (variableMap.count("n") != 1)
		{
			std::cout << "Please enter a single value for input `n'" << std::endl;
			return false;
		}
		out = variableMap["n"].as<std::size_t>();
		if (out <= 0)
		{
			std::cout << "Input `n' must be a positive integer" << std::endl;
			return false;
		}
		return true;
	}
	bool readInitialRadius(boost::program_options::variables_map& variableMap, int& out, std::string& message)
	{
		if(variableMap.count("initialRadius") != 1)
		{
			message = "Please enter a single value for input `initialRadius'";
			return false;
		}
		out = variableMap["initialRadius"].as<int>();
		if(out < 0)
		{
			message = "Input `initialRadius' must be a non-negative integer";
			return false;
		}
		return true;
	}

	bool readProbability(boost::program_options::variables_map& variableMap, double& out)
	{
		if(variableMap.count("opProbability") != 1)
		{
			std::cout << "Please enter a single value for input `opProbability'" << std::endl;
			return false;
		}
		double probability = variableMap["opProbability"].as<double>();
		if(probability > 1 || probability < 0)
		{
			std::cout << "Please enter a number between 0 and 1 for `opProbability'" << std::endl;
			return false;
		}
		out = probability;
		return true;
	}
	void readSeed(boost::program_options::variables_map& variableMap, boost::mt19937& randomSource)
	{
		if(variableMap.count("seed") > 0)
		{
			randomSource.seed(variableMap["seed"].as<int>());
		}
	}
}