#include <boost/program_options.hpp>
#include "computeConditionalProb.h"
#include "Arguments.h"
#include "includeMPFR.h"
#include "ArgumentsMPFR.h"
#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/algorithm/random_shuffle.hpp>
#include <algorithm>
#include "Context.h"
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/random_number_generator.hpp>
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
			("n", boost::program_options::value<int>(), "(int) The number of simulations to perform. ")
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
		int n;
		if(!readN(variableMap, n))
		{
			return 0;
		}
		mpfr_class opProbability;
		if(!readProbabilityString(variableMap, opProbability))
		{
			std::cout << "Error parsing numeric argument opProbability" << std::endl;
		}
		double opProbabilityD = opProbability.toDouble();
		Context context = Context::emptyContext();
		if(!readContext(variableMap, context, opProbabilityD))
		{
			return 0;
		}
		const std::vector<int>& interestVertices = context.getInterestVertices();

		boost::mt19937 randomSource;
		readSeed(variableMap, randomSource);

		const Context::internalGraph& graph = context.getGraph();
		const std::size_t nEdges = context.getNEdges();
		const std::size_t nVertices = boost::num_vertices(context.getGraph());

		std::vector<std::pair<int, int> > vertices;
		Context::internalGraph::edge_iterator current, end;
		boost::tie(current, end) = boost::edges(graph);
		for(;current != end; current++)
		{
			vertices.push_back(std::make_pair((int)current->m_source, (int)current->m_target));
		}

		std::vector<mpfr_class> exponentialRates(nEdges);
		{
			mpfr_class exponentialRate = -mpfr::log(1 - opProbability);
			std::fill(exponentialRates.begin(), exponentialRates.end(), exponentialRate);
		}
		//The sum of all the conditional probabilities (Used to get estimate and estimate of variance of estimate)
		mpfr_class sumConditional = 0, sumSquaredConditional = 0;
		//This stores the rates that go into the matrix exponential computation
		std::vector<mpfr_class> ratesForPMC;
		boost::random_number_generator<boost::mt19937> generator(randomSource);
		std::vector<int> edgeOrdering(boost::counting_iterator<int>(0), boost::counting_iterator<int>((int)nEdges));

		//The initial rate at the start
		mpfr_class sumAllRates = 0; 
		for (std::vector<mpfr_class>::iterator i = exponentialRates.begin(); i != exponentialRates.end(); i++)
		{
			sumAllRates += *i;
		}

		//The edges which are still under consideration
		std::vector<bool> alreadySeen(nEdges);
		//Only warn about stability once
		bool warnedStability = false;
		
		//scratch data used to compute conditional probabilities
		std::vector<mpfr_class> computeConditionalProbScratch;
		//Used to hold the connected component IDS
		std::vector<int> componentIDs(nVertices);
		for(int i = 0; i < n; i++)
		{
			//Simulate permutation
			boost::random_shuffle(edgeOrdering, generator);
			//No edges have yet been seen
			std::fill(alreadySeen.begin(), alreadySeen.end(), false);
			//The first rate is going to be this
			mpfr_class currentRate = sumAllRates;
			//which edge in the permutation are we currently looking at?
			int permutationCounter = 0;
			//these are going to be the rates for the matrix exponential
			ratesForPMC.clear();
			//reset connected component IDs
			componentIDs.clear();
			componentIDs.insert(componentIDs.begin(), boost::counting_iterator<int>(0), boost::counting_iterator<int>((int)nVertices));
			while(true)
			{
				//get out the edge index
				int edgeIndex = edgeOrdering[permutationCounter];
				//Has this edge been removed from consideration already?
				if(!alreadySeen[edgeIndex])
				{
					//If not, add the current rate
					ratesForPMC.push_back(currentRate);
					//Mark the edge as seen
					alreadySeen[edgeIndex] = true;
					//subtract the rate from the current rate
					currentRate -= exponentialRates[edgeIndex];
					//Get the old component IDs
					int firstComponentID = componentIDs[vertices[edgeIndex].first];
					int secondComponentID = componentIDs[vertices[edgeIndex].second];
					//go through all the edges that haven't been seen yet. 
					for(int j = 0; j < nEdges; j++)
					{
						if(!alreadySeen[j])
						{
							int firstVertex = vertices[j].first, secondVertex = vertices[j].second;
							if((componentIDs[firstVertex] == secondComponentID || componentIDs[firstVertex] == firstComponentID) && (componentIDs[secondVertex] == firstComponentID || componentIDs[secondVertex] == secondComponentID))
							{
								alreadySeen[j] = true;
								currentRate -= exponentialRates[j];
							}
						}
					}
					std::replace(componentIDs.begin(), componentIDs.end(), std::max(firstComponentID, secondComponentID), std::min(firstComponentID, secondComponentID));
					//determine whether or not we've hit the critical threshold
					bool connected = true;
					for(int k = 1; k < interestVertices.size(); k ++)
					{
						connected &= (componentIDs[interestVertices[k]] == componentIDs[interestVertices[0]]);
					}
					if(connected) break;
				}
				permutationCounter++;
			}
			mpfr_class additionalPart = computeConditionalProb(ratesForPMC, computeConditionalProbScratch);
			//mpfr_class additionalPart2 = computeConditionalProb(ratesForPMC);
			if(additionalPart > 1 && !warnedStability)
			{
				std::cout << "Numerical stability problem detected" << std::endl;
				warnedStability = true;
			}
			sumConditional += additionalPart;
			sumSquaredConditional += additionalPart*additionalPart;
		}
		mpfr_class estimateFirstMoment = sumConditional/n;
		mpfr_class estimateSecondMoment = sumSquaredConditional/n;
		mpfr_class varianceEstimate = estimateSecondMoment - estimateFirstMoment*estimateFirstMoment;
		mpfr_class sqrtVarianceEstimate = mpfr::sqrt(varianceEstimate/n);
		mpfr_class relativeErrorEstimate = sqrtVarianceEstimate / estimateFirstMoment;

		std::cout << "Unreliability probability estimate was " << estimateFirstMoment.toDouble() << std::endl;
		std::cout << "Relative error was " << relativeErrorEstimate.toDouble() << std::endl;
		return 0;
	}
}

int main(int argc, char **argv)
{
	return networkReliability::main(argc, argv);
}