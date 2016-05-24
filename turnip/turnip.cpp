#include <boost/program_options.hpp>
#include "computeConditionalProb.h"
#include "arguments.h"
#include "argumentsMPFR.h"
#include <boost/iterator/counting_iterator.hpp>
#include <algorithm>
#include "context.h"
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/random_number_generator.hpp>
#include "turnipImpl.h"
namespace networkReliability
{
	int main(int argc, char **argv)
	{
		mpfr_set_default_prec(1024);

		boost::program_options::options_description options("Usage");
		options.add_options()
			("completeGraph", boost::program_options::value<int>(), "(int) The number of vertices in the complete graph to use. Incompatible with graphFile and gridGraph. ")
			("gridGraph", boost::program_options::value<int>(), "(int) The dimension of the square grid graph to use. Incompatible with graphFile and completeGraph. ")
			("graphFile", boost::program_options::value<std::string>(), "(string) The path to a graphml file. Incompatible with gridGraph and completeGraph. ")
			("n", boost::program_options::value<std::size_t>(), "(int) The number of simulations to perform. ")
			("opProbability", boost::program_options::value<std::string>(), "(float) The probability that an edge is operational. ")
			("seed", boost::program_options::value<int>(), "(int) The random seed used to generate the random graphs. ")
			("interestVertices", boost::program_options::value<std::vector<int> >()->multitoken(), "(int) The vertices of interest, that should be connected. ")
			("help", "Display this message");
		
		boost::program_options::variables_map variableMap;
		try
		{
			boost::program_options::store(boost::program_options::parse_command_line(argc, argv, options), variableMap);
		}
		catch(boost::program_options::error& ee)
		{
			std::cerr << "Error parsing command line arguments: " << ee.what() << std::endl << std::endl;
			std::cerr << options << std::endl;
			return -1;
		}
		if(variableMap.count("help") > 0)
		{
			std::cout << 
				"This program estimates the probability that the given graph is unreliable using PMC. \n\n"
			;
			std::cout << options << std::endl;
			return 0;
		}
		std::size_t n;
		if(!readN(variableMap, n))
		{
			return 0;
		}
		mpfr_class opProbability;
		if(!readProbabilityString(variableMap, opProbability))
		{
			std::cout << "Error parsing numeric argument opProbability" << std::endl;
		}
		double opProbabilityD = opProbability.convert_to<double>();
		context contextObj = context::emptyContext();
		if(!readContext(variableMap, contextObj, opProbabilityD))
		{
			return 0;
		}

		boost::mt19937 randomSource;
		readSeed(variableMap, randomSource);

		const context::internalGraph& graph = contextObj.getGraph();
	
		turnipEqualRateInput input(randomSource, &graph, contextObj.getInterestVertices());
		//set up exponential rates
		input.exponentialRate = -boost::multiprecision::log(1 - opProbability);
		input.n = n;
		turnipEqualRate(input);
		if (input.warnedStability)
		{
			std::cout << "Numerical stability problem detected" << std::endl;
		}

		mpfr_class varianceEstimate = input.estimateSecondMoment - input.estimateFirstMoment*input.estimateFirstMoment;
		mpfr_class sqrtVarianceEstimate = boost::multiprecision::sqrt(varianceEstimate / input.n);
		mpfr_class relativeErrorEstimate = sqrtVarianceEstimate / input.estimateFirstMoment;

		std::cout << "Unreliability probability estimate was " << input.estimateFirstMoment.convert_to<double>() << std::endl;
		std::cout << "Relative error was " << relativeErrorEstimate.convert_to<double>() << std::endl;
		return 0;
	}
}

int main(int argc, char **argv)
{
	return networkReliability::main(argc, argv);
}
