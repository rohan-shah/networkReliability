#include "networkReliabilityObs.h"
#include <boost/program_options.hpp>
#include "arguments.h"
#include "argumentsMPFR.h"
#include "graphAlgorithms.h"
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/sum_kahan.hpp>
#include "formulaDriver.h"
#include "createFunctionBinary.h"
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
			("seed", boost::program_options::value<int>(), "(int) The random seed used to generate the random graphs. ")
			("interestVertices", boost::program_options::value<std::vector<int> >()->multitoken(), "(int) The vertices of interest, that should be connected. ")
			("function", boost::program_options::value<std::string>(), "The function to evaluate.")
			("functionFile", boost::program_options::value<std::string>(), "(string) File containing the function to evaluate")
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

		mpfr_class opProbability;
		if(!readProbabilityString(variableMap, opProbability))
		{
			std::cout << "Unable to read argument `opProbability'" << std::endl;
			return 0;
		}

		context contextObj = context::emptyContext();
		if(!readContext(variableMap, contextObj, opProbability))
		{
			return 0;
		}
		const std::size_t nEdges = contextObj.getNEdges();
		std::string function, message;
		formulaDriver driver(nEdges);
		if(!readFunctionFile(variableMap, driver, function, message))
		{
			std::cout << message << std::endl;
			return 0;
		}
		const std::size_t nFunctions = driver.result.size();
		boost::mt19937 randomSource;
		readSeed(variableMap, randomSource);

		boost::detail::depth_first_visit_restricted_impl_helper<context::internalGraph>::stackType stack;

		std::vector<int> components;
		std::vector<int> interestComponents;
		const std::vector<int> interestVertices = contextObj.getInterestVertices();
		std::vector<int> edgeValues(nEdges, 0);

		std::vector<boost::default_color_type> colorMap;
		std::vector<mpfr_class> sums(nFunctions, 0);
		long countDisconnected = 0;
		for(std::size_t i = 0; i < n; i++)
		{
			NetworkReliabilityObs obs(contextObj, randomSource);
			if(!isSingleComponent(contextObj, obs.getState(), components, stack, colorMap))
			{
				const edgeState* state = obs.getState();
				countDisconnected++;
				for(std::size_t functionCounter = 0; functionCounter < nFunctions; functionCounter++)
				{
					//For the edges which are relevant to the function, set those values in the vector to the state of the obseration
					for(std::vector<int>::iterator relevantEdgesIterator = driver.edgeIDs[functionCounter].begin(); relevantEdgesIterator != driver.edgeIDs[functionCounter].end(); relevantEdgesIterator++)
					{
						edgeValues[*relevantEdgesIterator] = ((state[*relevantEdgesIterator] & OP_MASK) > 0);
					}
					sums[functionCounter] += driver.result[functionCounter]->calculate(edgeValues);
				}
			}			
		}
		for(std::size_t functionCounter = 0; functionCounter < nFunctions; functionCounter++)
		{
			std::cout << "Estimated value of function " << functionCounter << " was " << (sums[functionCounter]/countDisconnected) << std::endl;
		}

		return 0;
	}
}
int main(int argc, char **argv)
{
	return networkReliability::main(argc, argv);
}
