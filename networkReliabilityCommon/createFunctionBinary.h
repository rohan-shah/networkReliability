#ifndef CREATE_FUNCTION_BINARY_HEADER_GUARD
#define CREATE_FUNCTION_BINARY_HEADER_GUARD
#include <string>
#include <vector>
namespace networkReliability
{
	bool createFunctionsBinary(std::string function, const std::vector<std::vector<int> >& edgeIDs, std::string& message, std::string& outputFile);
}
#endif
