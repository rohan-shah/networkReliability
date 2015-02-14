#include "Context.h"
#include <boost/program_options.hpp>
#include "Arguments.h"
#include "ArgumentsMPFR.h"
#include "EdgeState.h"
#include <boost/random/random_number_generator.hpp>
#include <boost/random/uniform_01.hpp>
#include "empiricalDistribution.h"
#include <boost/archive/binary_oarchive.hpp>
#include <fstream>
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
			("relativeError", boost::program_options::value<bool>()->default_value(false)->implicit_value(true), "(flag) Output relative error estimate")
			("standardDeviation", boost::program_options::value<bool>()->default_value(false)->implicit_value(true), "(flag) Output standard deviation estimate")
			("outputConditionalDistribution", boost::program_options::value<std::string>(), "(path) File to output the empirical conditional distribution")
			("optimiseMinCut", boost::program_options::value<bool>()->default_value(false)->implicit_value(true), "(flag) Should we try and optimise the number of mincut calls?")
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
		double inopProbabilityD = inopProbability.convert_to<double>();
		Context context = Context::emptyContext();
		if(!readContext(variableMap, context, opProbability))
		{
			return 0;
		}
		boost::mt19937 randomSource;
		readSeed(variableMap, randomSource);
		bool outputConditionalDistribution = variableMap.count("outputConditionalDistribution") > 0;
		bool relativeError = variableMap["relativeError"].as<bool>();
		bool standardDeviation = variableMap["standardDeviation"].as<bool>();
		const std::size_t nEdges = context.getNEdges();
		
		empiricalDistribution dist(true, nEdges, context);
		if(outputConditionalDistribution)
		{
			dist.hintDataCount(n);
		}
		bool optimiseMinCut = variableMap["optimiseMinCut"].as<bool>();
		if(context.getInterestVertices().size() > 2 && optimiseMinCut)
		{
			std::cout << "Can only specify option optimiseMinCut with 2-terminal reliability" << std::endl;
			return 0;
		}

		boost::random::uniform_01<float,float> uniformReal;


		//Vector used for mincut calculations
		std::vector<int> state(2*nEdges);
		//Similar vector of states, used to output the conditional distribution
		std::vector<EdgeState> outputStates(nEdges);
		//Sum over all the n simulations
		mpfr_class firstMoment = 0, secondMoment = 0;
		//likelihood ratio of current term
		mpfr_class currentLikelihoodRatio = 0;

		//Cache powers of the inopProbability 
		boost::scoped_array<mpfr_class> cachedInopPowers(new mpfr_class[nEdges+1]);
		for(std::size_t i = 0; i < nEdges+1; i++)
		{
			cachedInopPowers[i] = boost::multiprecision::pow(inopProbability, i);
		}
		//Get out the vector that holds the flow
		const std::vector<int>& residualCapacityVector = context.edgeResidualCapacityVector;
		//Cache the two mincut calculations that are the same every time - The first two. 
		mpfr_class minCutUpProb, minCutDownProb;
		std::fill(state.begin(), state.end(), 1);
		state[0] = state[1] = 0;
		int initialMinCutSizeDown = context.getMinCut(state);
		minCutDownProb = cachedInopPowers[initialMinCutSizeDown];

		state[0] = state[1] = HIGH_CAPACITY;
		int initialMinCutSizeUp = context.getMinCut(state);
		if(initialMinCutSizeUp >= HIGH_CAPACITY)
		{
			minCutUpProb = 0;
		}
		else minCutUpProb = cachedInopPowers[initialMinCutSizeUp];

		mpfr_class qTilde;
		mpfr_class initialQTilde = inopProbability*minCutDownProb;
		initialQTilde = initialQTilde / (initialQTilde  + opProbability * minCutUpProb);
		for(std::size_t i = 0; i < n; i++)
		{
			//we might break out of this loop early because the value of the indicator function is already known.
			bool exitedEarly = false;
			int indicatorValue;

			currentLikelihoodRatio = 0;
			std::fill(state.begin(), state.end(), 1);

			float random = uniformReal(randomSource);
			std::size_t edgeCounter = 0;
			if(random < initialQTilde)
			{
				state[0] = state[1] = 0; outputStates[0] = UNFIXED_INOP;
				currentLikelihoodRatio = inopProbability / initialQTilde;
				if(initialMinCutSizeDown == 0)
				{
					indicatorValue = 1;
					exitedEarly = true;
				}
			}
			else 
			{
				state[0] = state[1] = HIGH_CAPACITY; outputStates[0] = UNFIXED_OP;
				currentLikelihoodRatio = (1-inopProbability) / (1-initialQTilde);
				if(initialMinCutSizeUp >= HIGH_CAPACITY)
				{
					indicatorValue = 0;
					exitedEarly = true;
				}
			}
			if(!exitedEarly)
			{
				bool canReuseMinCut = false;
				int minCutSizeDown, minCutSizeUp;
				for(edgeCounter = 1; edgeCounter < nEdges; edgeCounter++)
				{
					state[2*edgeCounter] = state[2*edgeCounter+1] = 0;
					if(!canReuseMinCut || !optimiseMinCut) minCutSizeDown = context.getMinCut(state);
					minCutDownProb = cachedInopPowers[minCutSizeDown];

					state[2*edgeCounter] = state[2*edgeCounter+1] = HIGH_CAPACITY;
					if(!canReuseMinCut || !optimiseMinCut) minCutSizeUp = context.getMinCut(state);
					if(minCutSizeUp >= HIGH_CAPACITY)
					{
						minCutUpProb = 0;
					}
					else minCutUpProb = cachedInopPowers[ minCutSizeUp];

					qTilde = inopProbability * minCutDownProb;
					qTilde = qTilde / (qTilde  + opProbability * minCutUpProb);
					float random = uniformReal(randomSource);
					if(random < qTilde)
					{
						state[2*edgeCounter] = state[2*edgeCounter+1] = 0; outputStates[edgeCounter] = UNFIXED_INOP;
						currentLikelihoodRatio *= inopProbability / qTilde;
						if(minCutSizeDown == 0)
						{
							indicatorValue = 1;
							break;
						}
					}
					else 
					{
						state[2*edgeCounter] = state[2*edgeCounter+1] = HIGH_CAPACITY; outputStates[edgeCounter] = UNFIXED_OP;
						currentLikelihoodRatio *= opProbability / (1-qTilde);
						if(minCutSizeUp >= HIGH_CAPACITY)
						{
							indicatorValue = 0;
							break;
						}
					}
					canReuseMinCut = (residualCapacityVector[2*edgeCounter] == 1) && (residualCapacityVector[2*edgeCounter+1] == 1);
				}
			}
			firstMoment += indicatorValue * currentLikelihoodRatio;
			secondMoment += indicatorValue * currentLikelihoodRatio * currentLikelihoodRatio;
			//If we're outputting the conditional distribution then we need to simulate a state for the remaining edges
			if(outputConditionalDistribution && indicatorValue)
			{
				edgeCounter++;
				for(; edgeCounter < nEdges; edgeCounter++)
				{
					if(uniformReal(randomSource) < inopProbabilityD) outputStates[edgeCounter] = UNFIXED_INOP;
					else outputStates[edgeCounter] = UNFIXED_OP;
				}
				dist.add(&(outputStates[0]), currentLikelihoodRatio.convert_to<double>());
			}
		}
		mpfr_class estimate = (firstMoment / n);
		mpfr_class variance = (secondMoment/n) - estimate*estimate;
		mpfr_class sqrtVariance = boost::multiprecision::sqrt(variance / n);
		std::cout << "Estimated unreliability was " << estimate.convert_to<std::string>() << std::endl;
		if(relativeError)
		{
			mpfr_class relativeError = sqrtVariance / estimate;
			std::cout << "Estimated relative error was " << relativeError.convert_to<std::string>() << std::endl;
		}
		if(standardDeviation)
		{
			mpfr_class standardDeviation = boost::multiprecision::sqrt(variance);
			std::cout << "Standard deviation was " << standardDeviation << std::endl;
		}
		if(outputConditionalDistribution)
		{
			try
			{
				std::ofstream stream(variableMap["outputConditionalDistribution"].as<std::string>().c_str(), std::ios_base::binary);
				boost::archive::binary_oarchive archive(stream);
				archive << dist;
			}
			catch(std::runtime_error& err)
			{
				std::cout << "Error saving output distribution: " << err.what() << std::endl;
				return 0;
			}
		}
		return 0;
	}
}
int main(int argc, char** argv)
{
		return networkReliability::main(argc, argv);
}
