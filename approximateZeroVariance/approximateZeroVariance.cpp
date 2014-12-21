#include "Context.h"
#include <boost/program_options.hpp>
#include "Arguments.h"
#include "ArgumentsMPFR.h"
#include "EdgeState.h"
#include <boost/random/random_number_generator.hpp>
#include <boost/random/uniform_01.hpp>
namespace networkReliability
{
	int main(int argc, char** argv)
	{
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
				"This program estimates the probability that the given graph is unreliable using the approximate zero-variance importance sampling method of L'ecuyer (2011). \n\n"
			;
			std::cout << options << std::endl;
			return 0;
		}
		std::size_t n;
		if(!readN(variableMap, n))
		{
			return 0;
		}
		mpfr_class opProbability, inopProbability;
		if(!readProbabilityString(variableMap, opProbability))
		{
			std::cout << "Error parsing numeric argument opProbability" << std::endl;
			return 0;
		}
		inopProbability = 1 - opProbability;
		Context context = Context::emptyContext();
		if(!readContext(variableMap, context, opProbability))
		{
			return 0;
		}
		boost::mt19937 randomSource;
		readSeed(variableMap, randomSource);
		const std::size_t nEdges = context.getNEdges();
		
		boost::random::uniform_01<float,float> uniformReal;


		//Vector used for mincut calculations
		std::vector<int> state(2*nEdges);
		//Sum over all the n simulations
		mpfr_class sum = 0;
		//likelihood ratio of current term
		mpfr_class currentLikelihoodRatio = 0;

		//Cache the two mincut calculations that are the same every time - The first two. 
		mpfr_class minCutUpProb, minCutDownProb;
		std::fill(state.begin(), state.end(), 1);
		state[0] = state[1] = 0;
		int initialMinCutSizeDown = context.getMinCut(state);
		minCutDownProb = boost::multiprecision::pow(inopProbability, initialMinCutSizeDown);

		state[0] = state[1] = HIGH_CAPACITY;
		int initialMinCutSizeUp = context.getMinCut(state);
		if(initialMinCutSizeUp >= HIGH_CAPACITY)
		{
			minCutUpProb = 0;
		}
		else minCutUpProb = boost::multiprecision::pow(inopProbability, initialMinCutSizeUp);

		mpfr_class qTilde;
		mpfr_class initialQTilde = inopProbability*minCutDownProb;
		initialQTilde = initialQTilde / (initialQTilde  + opProbability * minCutUpProb);
		double tmp = initialQTilde.convert_to<double>();
		for(long i = 0; i < n; i++)
		{
			//we might break out of this loop early because the value of the indicator function is already known.
			bool exitedEarly = false;
			int indicatorValue;

			currentLikelihoodRatio = 0;
			std::fill(state.begin(), state.end(), 1);

			float random = uniformReal(randomSource);
			if(random < initialQTilde)
			{
				state[0] = state[1] = 0;
				currentLikelihoodRatio = inopProbability / initialQTilde;
				tmp = currentLikelihoodRatio.convert_to<double>();
				if(initialMinCutSizeDown == 0)
				{
					indicatorValue = 1;
					exitedEarly = true;
				}
			}
			else 
			{
				state[0] = state[1] = HIGH_CAPACITY;
				currentLikelihoodRatio = (1-inopProbability) / (1-initialQTilde);
				tmp = currentLikelihoodRatio.convert_to<double>();
				if(initialMinCutSizeUp >= HIGH_CAPACITY)
				{
					indicatorValue = 0;
					exitedEarly = true;
				}
			}
			if(!exitedEarly)
			{
				for(int edgeCounter = 1; edgeCounter < nEdges; edgeCounter++)
				{
					state[2*edgeCounter] = state[2*edgeCounter+1] = 0;
					int minCutSizeDown = context.getMinCut(state);
					minCutDownProb = boost::multiprecision::pow(inopProbability, minCutSizeDown);

					state[2*edgeCounter] = state[2*edgeCounter+1] = HIGH_CAPACITY;
					int minCutSizeUp = context.getMinCut(state);
					if(minCutSizeUp >= HIGH_CAPACITY)
					{
						minCutUpProb = 0;
					}
					else minCutUpProb = boost::multiprecision::pow(inopProbability, minCutSizeUp);

					qTilde = inopProbability * minCutDownProb;
					qTilde = qTilde / (qTilde  + opProbability * minCutUpProb);
					tmp = qTilde.convert_to<double>();
					float random = uniformReal(randomSource);
					if(random < qTilde)
					{
						state[2*edgeCounter] = state[2*edgeCounter+1] = 0;
						currentLikelihoodRatio *= inopProbability / qTilde;
						tmp = currentLikelihoodRatio.convert_to<double>();
						if(minCutSizeDown == 0)
						{
							indicatorValue = 1;
							break;
						}
					}
					else 
					{
						state[2*edgeCounter] = state[2*edgeCounter+1] = HIGH_CAPACITY;
						currentLikelihoodRatio *= opProbability / (1-qTilde);
						tmp = currentLikelihoodRatio.convert_to<double>();
						if(minCutSizeUp >= HIGH_CAPACITY)
						{
							indicatorValue = 0;
							break;
						}
					}
				}
			}
			sum += indicatorValue * currentLikelihoodRatio;
		}
		mpfr_class estimate = (sum / n);
		std::cout << "Estimated unreliability was " << estimate.convert_to<std::string>() << std::endl;
		return 0;
	}
}
int main(int argc, char** argv)
{
		return networkReliability::main(argc, argv);
}