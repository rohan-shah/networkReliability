#include <boost/lexical_cast.hpp>
#include "includeMPIRXX.h"
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include "Arguments.h"
#include "ArgumentsMPIR.h"
#include <boost/filesystem.hpp>
#include <locale>
#include <codecvt>
namespace networkReliability
{
	int main(int argc, char** argv)
	{
		mpf_set_default_prec(1024);

		boost::program_options::options_description exhaustiveProbabilityOptions("Usage");
		exhaustiveProbabilityOptions.add_options()
			("opProbability", boost::program_options::value<std::string>(), "(float) The probability that an edge is operational. ")
			("zeroProbability", boost::program_options::value<unsigned long long>(), "(Positive integer) Output the probability of getting a zero estimate from CMC, with specified number of samples");
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

		mpf_class probability_mpf, compProbability_mpf, one_mpf;
		mpf_set_d(one_mpf.get_mpf_t(), 1);
		if(!readProbabilityString(exhaustiveProbabilityVariableMap, probability_mpf))
		{
			return 0;
		}
		mpf_sub(compProbability_mpf.get_mpf_t(), one_mpf.get_mpf_t(), probability_mpf.get_mpf_t());

		std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> converter;
		std::vector<std::wstring> lines;
		std::string line;
		while(std::getline(std::cin, line))
		{
			lines.push_back(converter.from_bytes(line.c_str()));
		}
		std::size_t nLines = lines.size();
		if(nLines <= 2)
		{
			std::cout << "Input must contain more than two lines" << std::endl;
			return 0;
		}

		std::vector<std::wstring> splitFirstLine;
		boost::split(splitFirstLine, lines[0], boost::is_any_of("\""), boost::token_compress_on);
		if(splitFirstLine.size() != 5 || !boost::starts_with(splitFirstLine[0], L"Command ") || splitFirstLine[2] != L" with arguments " || splitFirstLine[4] != L"") 
		{
			std::cout << "First line was badly formatted" << std::endl;
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
		boost::scoped_array<mpf_class> sizeCounters(new mpf_class[nEdges+1]);
		for(int i = 0; i < nEdges+1; i++)
		{
			std::vector<std::string> splitCounterLine;
			boost::split(splitCounterLine, lines[i+2], boost::is_any_of(":"), boost::token_compress_on);
			boost::trim(splitCounterLine[0]);
			boost::trim(splitCounterLine[1]);
			if(splitCounterLine.size() != 2 || boost::lexical_cast<int>(splitCounterLine[0]) != i)
			{
				std::cout << "Invalid format for data line" << std::endl;
				return 0;
			}
			mpf_init_set_str(sizeCounters[i].get_mpf_t(), splitCounterLine[1].c_str(), 10);
		}

		mpf_class result = 0;
		for(int i = 0; i < nEdges+1; i++)
		{
			mpf_class probabilityPower;
			mpf_pow_ui(probabilityPower.get_mpf_t(), probability_mpf.get_mpf_t(), i);
			mpf_class compProbabilityPower;
			mpf_pow_ui(compProbabilityPower.get_mpf_t(), compProbability_mpf.get_mpf_t(), nEdges-i);

			result += sizeCounters[i] * compProbabilityPower * probabilityPower;
		}
		mpf_class unreliability;
		mpf_sub(unreliability.get_mpf_t(), one_mpf.get_mpf_t(), result.get_mpf_t());
		
		mp_exp_t exponent;
		char* resultCStr = mpf_get_str(NULL, &exponent, 10, 10, unreliability.get_mpf_t());
		std::string resultStr = resultCStr;
		std::cout << "Unreliability probability was " << resultStr.substr(0, 1) << "." << resultStr.substr(1, resultStr.size() - 1) << "e" << (exponent-1) << std::endl;
		free(resultCStr);

		if(exhaustiveProbabilityVariableMap.count("zeroProbability") == 1)
		{
			mpf_class reliability, reliabilityPower;
			mpf_sub(reliability.get_mpf_t(), one_mpf.get_mpf_t(), unreliability.get_mpf_t());
			mpf_pow_ui(reliabilityPower.get_mpf_t(), reliability.get_mpf_t(), exhaustiveProbabilityVariableMap["zeroProbability"].as<unsigned long long>());

			resultCStr = mpf_get_str(NULL, &exponent, 10, 10, reliabilityPower.get_mpf_t());
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