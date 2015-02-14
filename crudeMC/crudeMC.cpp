#include "NetworkReliabilityObs.h"
#include <boost/program_options.hpp>
#include "Arguments.h"
#include "ArgumentsMPFR.h"
#include "graphAlgorithms.h"
namespace networkReliability
{
	int main(int argc, char **argv)
	{
		boost::program_options::options_description options("Usage");
		options.add_options()
			("gridGraph", boost::program_options::value<int>(), "(int) The dimension of the square grid graph to use. Incompatible with graphFile and torusGraph. ")
			("graphFile", boost::program_options::value<std::string>(), "(string) The path to a graphml file. Incompatible with gridGraph")
			("completeGraph", boost::program_options::value<int>(), "(int) The number of vertices of the complete graph to use. ")
			("opProbability", boost::program_options::value<double>(), "(float) The probability that an edge is operational. ")
			("n", boost::program_options::value<std::size_t>(), "(int) The number of simulations to perform. ")
			("seed", boost::program_options::value<int>(), "(int) The random seed used to generate the random graphs. ")
			("interestVertices", boost::program_options::value<std::vector<int> >()->multitoken(), "(int) The vertices of interest, that should be connected. ")
			("splitting", boost::program_options::value<int>(), "(int) Should we estimating splitting level probabilities?")
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
				"This program estimates the probability that the given graph is unreliable for the given vertices. That is, if edges are retained with a certain probability, what is the probability that the specified vertices are not all in the same connected component?\n\n"
			;
			std::cout << options << std::endl;
			return 0;
		}

		std::size_t n;
		if(!readN(variableMap, n))
		{
			return 0;
		}

		double probability;
		if(!readProbability(variableMap, probability))
		{
			return 0;
		}

		Context context = Context::emptyContext();
		if(!readContext(variableMap, context, probability))
		{
			return 0;
		}

		boost::mt19937 randomSource;
		readSeed(variableMap, randomSource);

		boost::detail::depth_first_visit_restricted_impl_helper<Context::internalGraph>::stackType stack;

		std::vector<int> components;
		std::vector<int> interestComponents;
		const std::vector<int> interestVertices = context.getInterestVertices();

		std::vector<boost::default_color_type> colorMap;
		int countDisconnected = 0;

		for(std::size_t i = 0; i < n; i++)
		{
			NetworkReliabilityObs obs(context, randomSource);
			if(!isSingleComponent(context, obs.getState(), components, stack, colorMap))
			{
				countDisconnected++;
			}			
		}
		std::cout << "Estimated unreliability probability was " << countDisconnected << " / " << n << " = " << (float)countDisconnected / (float)n << std::endl;;

		return 0;
	}
}
int main(int argc, char **argv)
{
	return networkReliability::main(argc, argv);
}
