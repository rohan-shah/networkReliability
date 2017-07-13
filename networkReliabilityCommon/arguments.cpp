#include "arguments.h"
#include <boost/iterator/counting_iterator.hpp>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <iostream>
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
	bool readThresholds(boost::program_options::variables_map& variableMap, std::vector<double>& out, std::string& message)
	{
		if(variableMap.count("initialRadius") + variableMap.count("useSpatialDistances" ) != 1)
		{
			message = "Please enter exactly one of `initialRadius' and `useSpatialDistances'";
			return false;
		}
		else if(variableMap.count("initialRadius") == 1)
		{
			int initialRadius = variableMap["initialRadius"].as<int>();
			if(initialRadius < 0)
			{
				message = "Input `initialRadius' must be a non-negative integer";
				return false;
			}
			out.clear();
			for(int i = initialRadius; i > -1; i--) out.push_back(i);
			return true;
		}
		else
		{
			std::vector<double> thresholds = variableMap["useSpatialDistances"].as<std::vector<double> >();
			if(thresholds.size() != 2 && thresholds.size() != 3)
			{
				message = "Input spatial distances must consist of either two numbers or three; The two number case is a maximum distance and the number of steps to take. The three number case is a maximum distance, a minimum distance and the number of steps to take between the minimum and maximum (inclusive). Zero is added as the last threshold. ";
				return false;
			}
			if(thresholds.size() == 2)
			{
				double maximumDistance = thresholds[0], nSteps = thresholds[1];
				if(maximumDistance != maximumDistance || nSteps != nSteps || boost::math::isinf(nSteps) || boost::math::isinf(maximumDistance))
				{
					message = "Inputs for useSpatialDistances cannot be NA or infinity";
					return false;
				}
				if(abs(nSteps - floor(nSteps + 0.5)) > 1e-6)
				{
					message = "Number of steps to take for input `useSpatialDistances' must be an integer";
					return false;
				}
				if(nSteps < 0.5)
				{
					message = "Number of steps to take for input `useSpatialDistances' must be positive";
					return false;
				}
				if(nSteps == 1)
				{
					message = "Specifying useSpatialDistances with 1 step would be equivalent to crude Monte Carlo";
					return false;
				}
				int nStepsInt = (int)floor(nSteps + 0.5);
				out.clear();
				for(int currentStep = nStepsInt; currentStep > 0; currentStep--)
				{
					out.push_back(maximumDistance * ((double)(currentStep-1) / (double)(nStepsInt-1)));
				}
				return true;
			}
			else
			{
				double maximumDistance = thresholds[0], minimumDistance = thresholds[1], nSteps = thresholds[2];
				if(minimumDistance != minimumDistance || maximumDistance != maximumDistance || nSteps != nSteps || boost::math::isinf(nSteps) || boost::math::isinf(maximumDistance) || boost::math::isinf(minimumDistance))
				{
					message = "Inputs for useSpatialDistances cannot be NA or infinity";
					return false;
				}
				if(abs(nSteps - floor(nSteps + 0.5)) > 1e-6)
				{
					message = "Number of steps to take for input `useSpatialDistances' must be an integer";
					return false;
				}
				if(nSteps < 0.5)
				{
					message = "Number of steps to take for input `useSpatialDistances' must be positive";
					return false;
				}
				if(nSteps == 1)
				{
					message = "When specifying three inputs for useSpatialDistances, the number of steps must be at least 2";
					return false;
				}
				int nStepsInt = (int)floor(nSteps + 0.5);
				out.clear();
				for(int currentStep = nStepsInt; currentStep > 0; currentStep--)
				{
					out.push_back((maximumDistance - minimumDistance)* ((double)(currentStep-1) / (double)(nStepsInt-1)) + minimumDistance);
				}
				out.push_back(0);
				return true;
			}
		}
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
#ifdef HAS_BISON_AND_FLEX
	bool readFunctionFile(boost::program_options::variables_map& variableMap, formulaDriver& driver, std::string& function, std::string& message)
	{
		std::string functionFile;
		if(variableMap.count("function") + variableMap.count("functionFile") != 1)
		{
			message = "Exactly one of `function' or `functionFile' is requied";
			return false;
		}
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
					message = "Unable to read data from function file";
					return false;
			}
			function = std::string(std::istreambuf_iterator<char
>(functionFileStream), std::istreambuf_iterator<char>());
		}

		int parseResult = driver.parse(functionFile);
		//Unlike the temporary file, if we used one.
		if(variableMap.count("function")) unlink(functionFile.c_str());
		//Get out the unique edges that are of interest for this function.
		for(std::vector<std::vector<int> >::iterator i = driver.edgeIDs.begin(); i != driver.edgeIDs.end(); i++)
		{
			std::sort(i->begin(), i->end());
			i->erase(std::unique(i->begin(), i->end()), i->end());
		}
		if(parseResult != 0 || driver.message.size() != 0)
		{
		        message = "Errors parsing function: " + driver.message;
		        return false;
		}
		std::vector<std::string> functions;
		boost::split(functions, function, boost::is_any_of(","), boost::token_compress_on);
		if(functions.size() != driver.result.size() || functions.size() != driver.edgeIDs.size())
		{
			message = "Inconsistent number of functions";
			return false;
		}
		if(functions.size() == 0)
		{
			message = "Zero functions were input";
			return false;
		}
		return true;
	}
#endif
}
