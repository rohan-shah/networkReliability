#include <boost/graph/copy.hpp>
#include "NetworkReliabilityObs.h"
#include <boost/program_options.hpp>
#include "Arguments.h"
#include "ArgumentsMPFR.h"
#include "graphAlgorithms.h"
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/sum.hpp>
#include <boost/random/random_number_generator.hpp>
#include <boost/range/algorithm/random_shuffle.hpp>
#include "conditionalTurnip.h"
#include "conditionalPMC.h"
namespace networkReliability
{
	int main(int argc, char **argv)
	{
		boost::program_options::options_description options("Usage");
		options.add_options()
			("gridGraph", boost::program_options::value<int>(), "(int) The dimension of the square grid graph to use. Incompatible with graphFile and torusGraph. ")
			("graphFile", boost::program_options::value<std::string>(), "(string) The path to a graphml file. Incompatible with gridGraph")
			("completeGraph", boost::program_options::value<int>(), "(int) The number of vertices of the complete graph to use. ")
			("opProbability", boost::program_options::value<std::string>(), "(float) The probability that an edge is operational. ")
			("n", boost::program_options::value<std::size_t>(), "(int) The number of simulations to perform. ")
			("minimumFailedEdges", boost::program_options::value<int>(), "(int) Condition on having at least this number of edges failed")
			("seed", boost::program_options::value<int>(), "(int) The random seed used to generate the random graphs. ")
			("interestVertices", boost::program_options::value<std::vector<int> >()->multitoken(), "(int) The vertices of interest, that should be connected. ")
			("help", "Display this message");

		boost::program_options::variables_map variableMap;
		try
		{
			boost::program_options::store(boost::program_options::parse_command_line(argc, argv, options), variableMap);
		}
		catch (boost::program_options::error& ee)
		{
			std::cerr << "Error parsing command line arguments: " << ee.what() << std::endl << std::endl;
			std::cerr << options << std::endl;
			return -1;
		}
		if (variableMap.count("help") > 0)
		{
			std::cout <<
				"This program estimates the probability that the given graph is unreliable for the given vertices. That is, if edges are retained with a certain probability, what is the probability that the specified vertices are not all in the same connected component?\n\n"
				;
			std::cout << options << std::endl;
			return 0;
		}

		std::size_t n;
		if (!readN(variableMap, n))
		{
			return 0;
		}

		mpfr_class opProbability;
		if (!readProbabilityString(variableMap, opProbability))
		{
			return 0;
		}

		Context context = Context::emptyContext();
		if (!readContext(variableMap, context, opProbability))
		{
			return 0;
		}
		const std::size_t nEdges = context.getNEdges();

		boost::mt19937 randomSource;
		readSeed(variableMap, randomSource);

		if (variableMap.count("minimumFailedEdges") == 0)
		{
			std::cout << "Input minimumFailedEdges is required" << std::endl;
			return 0;
		}
		int minimumFailedEdges = variableMap["minimumFailedEdges"].as<int>();
		if (minimumFailedEdges < 0)
		{
			std::cout << "Input minimumFailedEdges must be non-negative" << std::endl;
			return 0;
		}

		const Context::internalGraph& originalGraph = context.getGraph();

		ConditionalTurnipInput input(randomSource, &originalGraph, context.getInterestVertices());
		input.minimumInoperative = minimumFailedEdges;
		input.exponentialRate = -boost::multiprecision::log(mpfr_class(1 - opProbability));
		input.n = n;
		input.edgeCounts.resize(nEdges, 1);

		conditionalPMC(input);
		mpfr_class estimatedVariance = input.estimateSecondMoment - input.estimateFirstMoment*input.estimateFirstMoment;
		std::cout << "Estimated conditional unreliability probability was " << input.estimateFirstMoment << std::endl;;
		std::cout << "Estimated relative error was " << boost::multiprecision::sqrt(estimatedVariance / n) / input.estimateFirstMoment << std::endl;
		return 0;
	}
}
int main(int argc, char **argv)
{
	return networkReliability::main(argc, argv);
}