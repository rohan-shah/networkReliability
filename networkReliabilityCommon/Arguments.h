#ifndef ARGUMENTS_HEADER_GUARD
#define ARGUMENTS_HEADER_GUARD
#include <boost/program_options.hpp>
#include <boost/random/mersenne_twister.hpp>
#include "Context.h"
namespace networkReliability
{
	bool readN(boost::program_options::variables_map& variableMap, std::size_t& out);
	bool readProbability(boost::program_options::variables_map& variableMap, double& out);
	void readSeed(boost::program_options::variables_map& variableMap, boost::mt19937& randomSource);
	bool readInitialRadius(boost::program_options::variables_map& variableMap, int& out, std::string& message);
	bool readFunctionFile(boost::program_options::variables_map& variableMap, std::string& functionFile, std::string& function, std::string& message);
}
#endif
