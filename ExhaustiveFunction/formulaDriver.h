#ifndef FORMULA_DRIVER_HEADER_GUARD
#define FORMULA_DRIVER_HEADER_GUARD
#include <string>
#include <map>
#include "formulaParser.hpp"
#include <istream>
#include <sstream>
#include <boost/shared_ptr.hpp>

#define YY_DECL \
  yy::formulaParser::symbol_type yylex(networkReliability::formulaDriver& driver)
YY_DECL;
namespace networkReliability
{
	class formulaDriver
	{	
	public:
		boost::shared_ptr<node> result;
		std::vector<int> edgeIDs;
		formulaDriver(int nEdges);
		virtual ~formulaDriver();
		std::string file;
		void scan_begin();
		void scan_end();
		bool trace_scanning;
		int parse(const std::string& data);
		bool trace_parsing;
		void error(const yy::location& l, const std::string& m);
		void error(const std::string& m);
		std::string message;
		int nEdges;
	};
}
#endif
