#ifndef ARGUMENTS_MPFR_HEADER_GUARD
#define ARGUMENTS_MPFR_HEADER_GUARD
#include <boost/program_options.hpp>
#include "includeMPFR.h"
#include "Context.h"
namespace networkReliability
{
	bool readProbabilityString(boost::program_options::variables_map& variableMap, mpfr_class& out);
	bool readContext(boost::program_options::variables_map& variableMap, Context& out, const mpfr_class& probability);
}
#endif