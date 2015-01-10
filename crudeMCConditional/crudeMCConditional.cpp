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

		mpfr_class probability;
		if (!readProbabilityString(variableMap, probability))
		{
			return 0;
		}

		Context context = Context::emptyContext();
		if (!readContext(variableMap, context, probability))
		{
			return 0;
		}
		const std::size_t nEdges = context.getNEdges();

		boost::mt19937 randomSource;
		readSeed(variableMap, randomSource);

		std::size_t minimumFailedEdges;
		if (variableMap.count("minimumFailedEdges") == 0)
		{
			minimumFailedEdges = context.getMinCutEdges();
		}
		else minimumFailedEdges = variableMap["minimumFailedEdges"].as<int>();
		if (minimumFailedEdges < 0)
		{
			std::cout << "Input minimumFailedEdges must be non-negative" << std::endl;
			return 0;
		}
		std::vector<int> permutation(boost::counting_iterator<int>(0), boost::counting_iterator<int>((int)nEdges));
		boost::random_number_generator<boost::mt19937> generator(randomSource);
		TruncatedBinomialDistribution::TruncatedBinomialDistribution operationalEdgesDist(nEdges, 0, nEdges - minimumFailedEdges, probability);
		boost::shared_array<EdgeState> edgeStates(new EdgeState[nEdges]);

		boost::detail::depth_first_visit_restricted_impl_helper<Context::internalGraph>::stackType stack;
		std::vector<int> components;
		std::vector<boost::default_color_type> colorMap;

		std::size_t countDisconnected = 0;
		for (std::size_t i = 0; i < n; i++)
		{
			boost::random_shuffle(permutation, generator);
			int nOperational = operationalEdgesDist(randomSource);
			for (int j = 0; j < nOperational; j++) edgeStates[permutation[j]] = UNFIXED_OP;
			for (int j = nOperational; j < nEdges; j++) edgeStates[permutation[j]] = UNFIXED_INOP;

			if (!isSingleComponent(context, edgeStates.get(), components, stack, colorMap))
			{
				countDisconnected++;
			}
		}
		std::cout << "Estimated conditional unreliability probability was " << countDisconnected << " / " << n << " = " << (float)countDisconnected / (float)n << std::endl;
		if (minimumFailedEdges <= context.getMinCutEdges())
		{
			TruncatedBinomialDistribution::TruncatedBinomialDistribution originalDistribution(nEdges, 0, nEdges, probability);
			mpfr_class conditioningProb = originalDistribution.getCumulativeProbability(nEdges - minimumFailedEdges);
			mpfr_class estimate = (conditioningProb * (float)countDisconnected / (float)n);
			mpfr_class variance = conditioningProb * conditioningProb * (float)countDisconnected * (float)(n - countDisconnected) / (float)(n*n*n);
			std::cout << "Estimated unreliability probability was " << estimate << std::endl;
			std::cout << "Relative error was " << (boost::multiprecision::sqrt(variance) / estimate) << std::endl;
		}
		return 0;
	}
}
int main(int argc, char **argv)
{
	return networkReliability::main(argc, argv);
}