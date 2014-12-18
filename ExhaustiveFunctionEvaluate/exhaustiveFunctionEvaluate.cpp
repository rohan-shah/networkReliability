#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "Arguments.h"
#include "ArgumentsMPFR.h"
#include <iomanip>
#include <fstream>
namespace networkReliability
{
	struct functionData
	{
		unsigned long long nEdges;
		long long functionValue;
		unsigned long long count;
	};
	int main(int argc, char **argv)
	{
		boost::program_options::options_description options("Usage");
		options.add_options()
			("opProbability", boost::program_options::value<std::string>(), "(float) The probability that an edge is operational") 
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
		mpfr_class opProbability;
		if(!readProbabilityString(variableMap, opProbability))
		{
			std::cout << "Unable to read input `opProbability'" << std::endl;
			return 0;
		}
		mpfr_class inopProbability = 1 - opProbability;
		std::vector<std::string> lines;
		std::string line;
		while(std::getline(std::cin, line))
		{
			lines.push_back(line);
		}
		const std::size_t nLines = lines.size();
		if(nLines <= 3)
		{
			std::cout << "Input must contain more than three lines" << std::endl;
			return 0;
		}
		std::vector<std::string> splitFirstLine;
		boost::split(splitFirstLine, lines[0], boost::is_any_of("\""), boost::token_compress_on);
		if(splitFirstLine.size() != 5 || !boost::starts_with(splitFirstLine[0], "Command ") || splitFirstLine[2] != " with arguments " || splitFirstLine[4] != "")
		{
			std::cout << "First line was badly formatted" << std::endl;
			return 0;
		}

		//Now look at second line
		int nFunctions;
		int result = sscanf(lines[1].c_str(), "%d functions input", &nFunctions);
		if(result != 1)
		{
			std::cout << "Second line was badly formatted" << std::endl;
		}

		//Now look at third line
		int nEdges;
		result = sscanf(lines[2].c_str(), "Graph had %d edges", &nEdges);
		if(result != 1)
		{
			std::cout << "Third line was badly formatted" << std::endl;
			return 0;
		}
		int lineCounter = 3;
		for(int functionCounter = 0; functionCounter < nFunctions; functionCounter++)
		{
			int nValues;
			result = sscanf(lines[lineCounter++].c_str(), ("Function " + boost::lexical_cast<std::string>(functionCounter) + " took on %d values").c_str(), &nValues);
			if(result != 1)
			{
				std::cout << "Line " << lineCounter << " was badly formatted" << std::endl;
				return 0;
			}
			std::vector<functionData> allFunctionData;
			std::vector<std::string>::iterator endIterator = lines.begin() + lineCounter + nValues * nEdges;
			for(std::vector<std::string>::iterator i = lines.begin() + lineCounter; i != endIterator; i++)
			{
				functionData readData;
				result = sscanf(i->c_str(), "%llu edges, value %lli : %llu", &readData.nEdges, &readData.functionValue, &readData.count);
				if(result != 3)
				{
					std::cout << "Input data line was badly formatted" << std::endl;
					return 0;
				}
				allFunctionData.push_back(readData);
			}
			mpfr_class expectedValue = 0, probabilitySum = 0;
			for(std::vector<functionData>::iterator i = allFunctionData.begin(); i != allFunctionData.end(); i++)
			{
				expectedValue += boost::multiprecision::pow(opProbability, i->nEdges) * boost::multiprecision::pow(inopProbability, nEdges - i->nEdges) * i->count * i->functionValue;
				probabilitySum += boost::multiprecision::pow(opProbability, i->nEdges) * boost::multiprecision::pow(inopProbability, nEdges - i->nEdges) * i->count;
			}
			mpfr_class scaledExpectedValue = expectedValue/probabilitySum;
			std::cout << "Expected value of function " << functionCounter << " was " << scaledExpectedValue.str() << std::endl;
			lineCounter += nValues*(nEdges+1);
		}
		return 0;
	}
}
int main(int argc, char **argv)
{
	return networkReliability::main(argc, argv);
}
