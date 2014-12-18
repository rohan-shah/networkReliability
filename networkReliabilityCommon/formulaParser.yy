%skeleton "lalr1.cc"
%defines
%define parser_class_name { formulaParser }
%define api.token.constructor
%define api.value.type variant
%define parse.assert
%code requires
{
	#include <string.h>
	#define YYDEBUG 1
	#include "astCode.hpp"
	namespace networkReliability
	{
		class formulaDriver;
	}
}
%param { networkReliability::formulaDriver& driver }
%locations
%initial-action
{
	@$.begin.filename = @$.end.filename = &driver.file;
};
%define parse.trace
%define parse.error verbose
%printer { std::cout << $$ << " "; } <int>;
%verbose
%define api.token.prefix {TOK_}
%code
{
	#include "formulaDriver.h"
}
%token 
	END 0 	"end of file"
	MINUS	"-"
	PLUS	"+"
	DIVIDE	"/"
	TIMES	"*"
	LPAREN	"("
	RPAREN	")"
	COMMA	","
	INTEGER
	EDGE
;
%type <node*> exp
%type <int> INTEGER EDGE
%left MINUS PLUS
%left TIMES DIVIDE
%nonassoc NEG
%nonassoc LPAREN RPARE

%%
%start unit;
unit: 	exp
		{ 
			driver.result.push_back(boost::shared_ptr<node>($1));
			driver.edgeIDs.push_back(std::move(driver.currentEdgeIDs));
			driver.currentEdgeIDs.clear();
		}
	| unit COMMA exp 		
		{
			driver.result.push_back(boost::shared_ptr<node>($3)); 
			driver.edgeIDs.push_back(std::move(driver.currentEdgeIDs));
			driver.currentEdgeIDs.clear();
		}
exp:
	  exp PLUS exp		{ $$ = new plusNode($1, $3); }
	| exp MINUS exp		{ $$ = new minusNode($1, $3); }
	| exp TIMES exp		{ $$ = new timesNode($1, $3); }
	| exp DIVIDE exp	{ $$ = new divideNode($1, $3); }
	| MINUS exp %prec NEG	{ $$ = new negateNode($2); }
	| LPAREN exp RPAREN	{ $$ = $2; };
	| INTEGER		{ $$ = new integerNode($1); }
	| EDGE
		{
			if($1 >= driver.nEdges || $1 < 0)
			{
				std::stringstream ss;
				ss << "Edge " << $1 << " is out of bounds";
				error(@$, ss.str());
				$$ = NULL;
			}
			driver.currentEdgeIDs.push_back($1);
			$$ = new edgeNode($1);
		}
%%
	void yy::formulaParser::error(const location_type& l, const std::string& m)
	{
		driver.error(l,m);
	}
