#include "ArgumentsMPIR.h"
namespace networkReliability
{
	bool readProbabilityString(boost::program_options::variables_map& variableMap, mpf_class& out)
	{
		int retVal = mpf_set_str(out.get_mpf_t(), variableMap["opProbability"].as<std::string>().c_str(), 10);
		return retVal == 0;
	}
}