#ifndef ARGUMENTS_MPFR_HEADER_GUARD
#define ARGUMENTS_MPFR_HEADER_GUARD
#include <boost/program_options.hpp>
#include "includeMPFR.h"
namespace networkReliability
{
	bool readProbabilityString(boost::program_options::variables_map& variableMap, mpfr_class& out);
}
#endif