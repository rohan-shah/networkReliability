#include <boost/lexical_cast.hpp>
#include "includeMPFR.h"
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include "Arguments.h"
#include "ArgumentsMPFR.h"
#include <boost/filesystem.hpp>
#include <locale>
#include <boost/math/distributions.hpp>
namespace networkReliability
{
	int main(int argc, char** argv)
	{
		mpfr_set_default_prec(1024);

		boost::program_options::options_description exhaustiveProbabilityOptions("Usage");
		exhaustiveProbabilityOptions.add_options()
			("opProbability", boost::program_options::value<std::string>(), "(float) The probability that an edge is operational. ")
			("minimumFailedEdges", boost::program_options::value<int>(), "(int) Condition on having at least this number of edges failed")
			("zeroProbability", boost::program_options::value<unsigned long long>(), "(Positive integer) Output the probability of getting a zero estimate from CMC, with specified number of samples")
			("help", "Display this message");
		boost::program_options::variables_map exhaustiveProbabilityVariableMap;

		try
		{
			boost::program_options::store(boost::program_options::parse_command_line(argc, argv, exhaustiveProbabilityOptions), exhaustiveProbabilityVariableMap);
		}
		catch(boost::program_options::error& ee)
		{
			std::cerr << "Error parsing command line arguments from input file: " << ee.what() << std::endl << std::endl;
			std::cerr << exhaustiveProbabilityOptions << std::endl;
			return -1;
		}
		if(exhaustiveProbabilityVariableMap.count("help") > 0)
		{
			std::cout << 
				"This program computes the probability that the given graph is unreliable for the given vertices. That is, if edges are retained with a certain probability, what is the probability that the specified vertices are not all in the same connected component? Standard input must be redirected from a file created by ExhaustiveSearch.exe" << std::endl << std::endl;
			;
			std::cout << exhaustiveProbabilityOptions << std::endl;
			return 0;
		}

		if(exhaustiveProbabilityVariableMap.count("zeroProbability") > 1)
		{
			std::cout << "Multiple values found for input `zeroProbability'" << std::endl;
			return 0;
		}

		mpfr_class probability_mpfr;
		if(!readProbabilityString(exhaustiveProbabilityVariableMap, probability_mpfr))
		{
			return 0;
		}
		mpfr_class compProbability_mpfr = 1 - probability_mpfr;

		std::vector<std::string> lines;
		std::string line;
		while(std::getline(std::cin, line))
		{
			boost::algorithm::trim(line);
			lines.push_back(line);
		}
		std::size_t nLines = lines.size();
		if(nLines <= 2)
		{
			std::cout << "Input must contain more than two lines" << std::endl;
			return 0;
		}

		std::vector<std::string> splitFirstLine;
		boost::split(splitFirstLine, lines[0], boost::is_any_of("\""), boost::token_compress_on);
		if(splitFirstLine.size() != 5 || !boost::starts_with(splitFirstLine[0], "Command ") || splitFirstLine[2] != " with arguments " || splitFirstLine[4] != "") 
		{
			std::cout << "First line was badly formatted" << std::endl;
			return 0;
		}

		if (lines[1] != "Number of connected subgraphs with that number of edges")
		{
			std::cout << "Second line was badly formatted" << std::endl;
			return 0;
		}

		/*boost::system::error_code ec;
		boost::filesystem::path workingDirectory(splitFirstLine[1]);
		boost::filesystem::current_path(workingDirectory, ec);
		if(ec)
		{
			std::cout << "Error setting directory" << std::endl;
			return 0;
		}

		boost::program_options::options_description exhaustiveSearchOptions("Usage");
		exhaustiveSearchOptions.add_options()
			("gridGraph", boost::program_options::value<int>(), "(int) The dimension of the square grid graph to use. Incompatible with graphFile and torusGraph. ")
			("graphFile", boost::program_options::value<std::string>(), "(string) The path to a graphml file. Incompatible with gridGraph")
			("completeGraph", boost::program_options::value<int>(), "(int) The number of vertices of the complete graph to use. ")
			("interestVertices", boost::program_options::value<std::vector<int> >()->multitoken(), "(int) The vertices of interest, that should be connected. ");

		boost::program_options::variables_map exhaustiveSearchVariableMap;
		try
		{
			std::vector<std::wstring> argumentsVector;
			boost::char_separator<wchar_t> sep(L" \n\r");
			boost::tokenizer<boost::char_separator<wchar_t>, std::wstring::const_iterator, std::wstring> tok(splitFirstLine[3], sep);
			std::copy(tok.begin(), tok.end(), std::back_inserter(argumentsVector));
			boost::program_options::store(boost::program_options::basic_command_line_parser<wchar_t>(argumentsVector).options(exhaustiveSearchOptions).run(), exhaustiveSearchVariableMap);
		}
		catch(boost::program_options::error& ee)
		{
			std::cerr << "Error parsing command line arguments from input file: " << ee.what() << std::endl << std::endl;
			std::cerr << exhaustiveSearchOptions << std::endl;
			return -1;
		}

		Context context = Context::emptyContext();
		std::string message;
		if(!readContext(exhaustiveSearchVariableMap, context, 0.5))
		{
			std::cout << "Internal error, could not read original graph: " << message << std::endl;
			return 0;
		}

		const std::size_t nEdges = context.getNEdges();
		if(lines.size() != nEdges+3)
		{
			std::cout << "Wrong number of lines in input file" << std::endl;
			return 0;
		}*/
		const std::size_t nEdges = lines.size() - 3;

		typedef long long counterType;
		boost::scoped_array<mpfr_class> sizeCounters(new mpfr_class[nEdges+1]);
		for(std::size_t i = 0; i < nEdges+1; i++)
		{
			std::vector<std::string> splitCounterLine;
			boost::split(splitCounterLine, lines[i+2], boost::is_any_of(":"), boost::token_compress_on);
			boost::trim(splitCounterLine[0]);
			boost::trim(splitCounterLine[1]);
			if(splitCounterLine.size() != 2 || boost::lexical_cast<std::size_t>(splitCounterLine[0]) != i)
			{
				std::cout << "Invalid format for data line" << std::endl;
				return 0;
			}
			sizeCounters[i] = mpfr_class(splitCounterLine[1]);
		}

		int minimumFailedEdges;
		if (exhaustiveProbabilityVariableMap.count("minimumFailedEdges") > 0)
		{
			minimumFailedEdges = exhaustiveProbabilityVariableMap["minimumFailedEdges"].as<int>();
			if (minimumFailedEdges < 0)
			{
				std::cout << "Input minimumFailedEdges must be non-negative" << std::endl;
				return 0;
			}
		}
		else minimumFailedEdges = 0;

		mpfr_class result = 0;
		for(std::size_t i = 0; i < nEdges+1-minimumFailedEdges; i++)
		{
			mpfr_class probabilityPower = boost::multiprecision::pow(probability_mpfr, i);
			mpfr_class compProbabilityPower = boost::multiprecision::pow(compProbability_mpfr, (unsigned long)(nEdges - i));

			result += sizeCounters[i] * compProbabilityPower * probabilityPower;
		}

		if (minimumFailedEdges > 0)
		{
			boost::math::binomial_distribution<mpfr_class> dist(nEdges, probability_mpfr);
			result /= boost::math::cdf(dist, nEdges - minimumFailedEdges);
		}
		mpfr_class unreliability = 1 - result;
		
		mp_exp_t exponent;
		char* resultCStr = mpfr_get_str(NULL, &exponent, 10, 10, unreliability.backend().data(), MPFR_RNDN);
		std::string resultStr = resultCStr;
		std::cout << "Unreliability probability was " << resultStr.substr(0, 1) << "." << resultStr.substr(1, resultStr.size() - 1) << "e" << (exponent-1) << std::endl;
		free(resultCStr);

		if(exhaustiveProbabilityVariableMap.count("zeroProbability") == 1)
		{
			mpfr_class reliability = 1 - unreliability;
			mpfr_class reliabilityPower = boost::multiprecision::pow(reliability, (unsigned long)exhaustiveProbabilityVariableMap["zeroProbability"].as<unsigned long long>());

			resultCStr = mpfr_get_str(NULL, &exponent, 10, 10, reliabilityPower.backend().data(), MPFR_RNDN);
			resultStr = resultCStr;
			std::cout << "Probability of estimating 0 for " << exhaustiveProbabilityVariableMap["zeroProbability"].as<unsigned long long>() << " CMC samples is " << resultStr.substr(0, 1) << "." << resultStr.substr(1, resultStr.size() - 1) << "e" << (exponent - 1) << std::endl;
			free(resultCStr);
		}
		return 0;
	}
}
int main(int argc, char** argv)
{
	return networkReliability::main(argc, argv);

}
