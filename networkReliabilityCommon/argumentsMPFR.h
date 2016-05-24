#ifndef ARGUMENTS_MPFR_HEADER_GUARD
#define ARGUMENTS_MPFR_HEADER_GUARD
#include <boost/program_options.hpp>
#include "includeMPFRNetworkReliability.h"
#include "context.h"
namespace networkReliability
{
	bool readProbabilityString(boost::program_options::variables_map& variableMap, mpfr_class& out);
	bool readContext(boost::program_options::variables_map& variableMap, context& out, const mpfr_class& probability);
}
#endif
