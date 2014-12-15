%{
#include <string.h>
#include "formulaParser.hpp"
#include "formulaDriver.h"
static yy::location loc;
%}
%option noyywrap nounput batch debug noinput

%{
#define YY_USER_ACTION loc.columns(yyleng);
%}

blank	[ \t]
integer	[0-9]+
edge	e[0-9]+
%%
%{
loc.step();
%}

{blank}+ 	loc.step();
[\n]+		loc.lines(yyleng); loc.step();
"-"		return yy::formulaParser::make_MINUS(loc);
"+"		return yy::formulaParser::make_PLUS(loc);
"/"		return yy::formulaParser::make_DIVIDE(loc);
"*"		return yy::formulaParser::make_TIMES(loc);
"("		return yy::formulaParser::make_LPAREN(loc);
")"		return yy::formulaParser::make_RPAREN(loc);
{integer}	{
			errno = 0;
			long n = strtol(yytext, NULL, 10);
			if(errno != 0) driver.error(loc, "Unable to parse integer");
			return yy::formulaParser::make_INTEGER((int)n, loc);
		}
{edge}		{
			errno = 0;
			long n = strtol(yytext+1, NULL, 10);
			if(errno != 0) driver.error(loc, "Unable to parse edge number");
			return yy::formulaParser::make_EDGE((int)n, loc);
		}
.		driver.error(loc, "invalid character");
<<EOF>>		return yy::formulaParser::make_END(loc);
%%
namespace networkReliability
{
	void formulaDriver::scan_begin()
	{
		yy_flex_debug = trace_scanning;
		if(!(yyin = fopen(file.c_str(), "r")))
		{
			message += "Unable to open file " + file;
		}
	}
	void formulaDriver::scan_end()
	{
		fclose(yyin);
	}
}
