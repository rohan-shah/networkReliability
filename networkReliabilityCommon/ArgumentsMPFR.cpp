#include "ArgumentsMPFR.h"
namespace networkReliability
{
	bool readProbabilityString(boost::program_options::variables_map& variableMap, mpfr_class& out)
	{
		if(variableMap.count("opProbability") < 1) return false;
		out = mpfr_class(variableMap["opProbability"].as<std::string>());
		return true;
	}
	bool readContext(boost::program_options::variables_map& variableMap, Context& out, const mpfr_class& probability)
	{
		boost::shared_ptr<std::vector<int> > interestVertices;
		if (variableMap.count("interestVertices") != 1)
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
		if (minInterest < 0)
		{
			std::cout << "Input `interestVertices' cannot contain negative indices" << std::endl;
			return false;
		}
		if (variableMap.count("graphFile") + variableMap.count("gridGraph") + variableMap.count("torusGraph") + variableMap.count("completeGraph") != 1)
		{
			std::cout << "Please enter exactly one of `completeGrahp', `gridGraph', `graphFile' or `torusGraph'" << std::endl;
			return false;
		}
		else if (variableMap.count("graphFile") == 1)
		{
			if (interestVertices->size() <= 1)
			{
				std::cout << "Input `interestVertices' must contain at least two vertex indices" << std::endl;
				return false;
			}

			bool successful;
			std::string message;
			out = Context::fromFile(variableMap["graphFile"].as<std::string>(), successful, interestVertices, message, probability);
			if (!successful)
			{
				std::cout << "Error reading graphml file. " << message << ". Exiting..." << std::endl;
				return false;
			}
			std::size_t nVertices = boost::num_vertices(out.getGraph());
		}
		else if (variableMap.count("gridGraph") == 1)
		{
			if (interestVertices->size() <= 1)
			{
				std::cout << "Input `interestVertices' must contain at least two vertex indices" << std::endl;
				return false;
			}

			int gridDimension = variableMap["gridGraph"].as<int>();
			if (gridDimension <= 0)
			{
				std::cout << "Input `gridGraph' must be a positive integer" << std::endl;
				return false;
			}
			if (maxInterest >= gridDimension*gridDimension || minInterest < 0)
			{
				std::cout << "Input 'interestVertices' must contain numbers between 0 and (nVertices - 1) inclusive" << std::endl;
				return false;
			}
			out = Context::gridContext(gridDimension, interestVertices, probability);
		}
		else if (variableMap.count("torusGraph") == 1)
		{
			std::cout << "Input `torusGraph' not yet supported" << std::endl;
			return false;
		}
		else if (variableMap.count("completeGraph") == 1)
		{
			int nVertices = variableMap["completeGraph"].as<int>();
			if (interestVertices->size() != 1)
			{
				std::cout << "For a complete graph, input interestVertices must contain the number of vertices of interest" << std::endl;
				return false;
			}
			if (minInterest < 2)
			{
				std::cout << "There must be at least two vertices of interest" << std::endl;
				return false;
			}
			if (minInterest > nVertices)
			{
				std::cout << "Input interestVertices was too large" << std::endl;
				return false;
			}
			out = Context::completeContext(nVertices, minInterest, probability);
		}
		if(variableMap.count("useMinCut") > 0) out.setMinCut(variableMap["useMinCut"].as<bool>());
		return true;
	}

}
