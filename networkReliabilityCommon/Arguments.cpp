#include "Arguments.h"
#include <boost/iterator/counting_iterator.hpp>
namespace networkReliability
{
	bool readN(boost::program_options::variables_map& variableMap, int& out)
	{
		if(variableMap.count("n") != 1)
		{
			std::cout << "Please enter a single value for input `n'" << std::endl;
			return false;
		}
		out = variableMap["n"].as<int>();
		if(out <= 0)
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
	bool readContext(boost::program_options::variables_map& variableMap, Context& out, double probability)
	{
		boost::shared_ptr<std::vector<int> > interestVertices;
		if(variableMap.count("interestVertices") != 1)
		{
			std::cout << "Please enter only one value for input `interestVertices'" << std::endl;
			return false;
		}
		{
			std::vector<int> tmp = variableMap["interestVertices"].as<std::vector<int> >();
			interestVertices = boost::shared_ptr<std::vector<int> >(new std::vector<int>());
			interestVertices->insert(interestVertices->begin(), tmp.begin(), tmp.end());
		}

		int minInterest = *std::min_element(interestVertices->begin(), interestVertices->end());
		int maxInterest = *std::max_element(interestVertices->begin(), interestVertices->end());
		if(minInterest < 0)
		{
			std::cout << "Input `interestVertices' cannot contain negative indices" << std::endl;
			return false;
		}
		if(variableMap.count("graphFile") + variableMap.count("gridGraph") + variableMap.count("torusGraph") + variableMap.count("completeGraph") != 1)
		{
			std::cout << "Please enter exactly one of `completeGrahp', `gridGraph', `graphFile' or `torusGraph'" << std::endl;
			return false;
		}
		else if(variableMap.count("graphFile") == 1)
		{
			if(interestVertices->size() <= 1)
			{
				std::cout << "Input `interestVertices' must contain at least two vertex indices" << std::endl;
				return false;
			}

			bool successful;
			std::string message;
			out = Context::fromFile(variableMap["graphFile"].as<std::string>(), successful, interestVertices, message, probability);
			if(!successful)
			{
				std::cout << "Error reading graphml file. " << message << ". Exiting..." << std::endl;
				return false;
			}
			std::size_t nVertices = boost::num_vertices(out.getGraph());
		}
		else if(variableMap.count("gridGraph") == 1)
		{
			if(interestVertices->size() <= 1)
			{
				std::cout << "Input `interestVertices' must contain at least two vertex indices" << std::endl;
				return false;
			}

			int gridDimension = variableMap["gridGraph"].as<int>();
			if(gridDimension <= 0)
			{
				std::cout << "Input `gridGraph' must be a positive integer" << std::endl;
				return false;
			}
			if(maxInterest >= gridDimension*gridDimension || minInterest < 0)
			{
				std::cout << "Input 'interestVertices' must contain numbers between 0 and (nVertices - 1) inclusive" << std::endl;
				return false;
			}
			out = Context::gridContext(gridDimension, interestVertices, probability);
		}
		else if(variableMap.count("torusGraph") == 1)
		{
			std::cout << "Input `torusGraph' not yet supported" << std::endl;
			return false;
		}
		else if(variableMap.count("completeGraph") == 1)
		{
			int nVertices = variableMap["completeGraph"].as<int>();
			if(interestVertices->size() != 1)
			{
				std::cout << "For a complete graph, input interestVertices must contain the number of vertices of interest" << std::endl;
				return false;
			}
			if(minInterest < 2)
			{
				std::cout << "There must be at least two vertices of interest" << std::endl;
				return false;
			}
			if(minInterest > nVertices)
			{
				std::cout << "Input interestVertices was too large" << std::endl;
				return false;
			}
			out = Context::completeContext(nVertices, minInterest, probability);
		}
		if(variableMap.count("distributionsCache") > 0)
		{
			std::vector<std::string> distributionsCache = variableMap["distributionsCache"].as<std::vector<std::string> >();
			for(std::vector<std::string>::const_iterator i = distributionsCache.begin(); i != distributionsCache.end(); i++)
			{
				out.loadDistributions(*i);
			}
		}
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