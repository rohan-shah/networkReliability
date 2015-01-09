#ifndef ARGUMENTS_HEADER_GUARD
#define ARGUMENTS_HEADER_GUARD
#include <boost/program_options.hpp>
#include <boost/random/mersenne_twister.hpp>
#include "Context.h"
#ifdef HAS_BISON_AND_FLEX
	#include "formulaDriver.h"
#endif
namespace networkReliability
{
	bool readN(boost::program_options::variables_map& variableMap, std::size_t& out);
	bool readProbability(boost::program_options::variables_map& variableMap, double& out);
	void readSeed(boost::program_options::variables_map& variableMap, boost::mt19937& randomSource);
	bool readThresholds(boost::program_options::variables_map& variableMap, std::vector<double>& out, std::string& message);
#ifdef HAS_BISON_AND_FLEX
	bool readFunctionFile(boost::program_options::variables_map& variableMap, formulaDriver& driver, std::string& function, std::string& message);
#endif
}
#endif
