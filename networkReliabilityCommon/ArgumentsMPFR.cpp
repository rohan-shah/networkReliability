#include "ArgumentsMPFR.h"
namespace networkReliability
{
	bool readProbabilityString(boost::program_options::variables_map& variableMap, mpfr_class& out)
	{
		out = mpfr_class(variableMap["opProbability"].as<std::string>());
		return true;
	}
}