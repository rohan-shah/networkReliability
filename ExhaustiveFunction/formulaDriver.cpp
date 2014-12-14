#include "formulaDriver.h"
#include "formulaParser.hpp"
#include <sstream>
namespace networkReliability
{
	formulaDriver::formulaDriver(int nEdges)
		:trace_scanning(false), trace_parsing(false), message(""), nEdges(nEdges)
	{
	}
	formulaDriver::~formulaDriver()
	{
	}
	int formulaDriver::parse(const std::string& file)
	{
		this->file = file;
		message = "";
		scan_begin();
		if(message != "") return 1;

		yy::formulaParser parser(*this);
		parser.set_debug_level(trace_parsing);
		int res  = parser.parse();
		scan_end();
		return res;
	}
	void formulaDriver::error(const yy::location& l, const std::string& m)
	{
		std::stringstream ss;
		ss << l << ": " << m << std::endl;;
		message += ss.str();
	}
	void formulaDriver::error(const std::string& message)
	{
		std::stringstream ss;
		ss << message << std::endl;
		this->message += ss.str();
	}
}
