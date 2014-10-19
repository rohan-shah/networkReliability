#include "ArgumentsMPFR.h"
namespace networkReliability
{
	bool readProbabilityString(boost::program_options::variables_map& variableMap, mpfr_class& out)
	{
		int retVal = mpfr_set_str(out.mpfr_ptr(), variableMap["opProbability"].as<std::string>().c_str(), 10, MPFR_RNDN);
		return retVal == 0;
	}
}