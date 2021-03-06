#include "conditionalPMC.h"
#include <boost/iterator/counting_iterator.hpp>
#include <boost/random/random_number_generator.hpp>
#include <boost/range/algorithm/random_shuffle.hpp>
#include <boost/random/exponential_distribution.hpp>
#include "computeConditionalProb.h"
namespace networkReliability
{
	ConditionalTurnipInput::ConditionalTurnipInput(boost::mt19937& randomSource, const context::internalGraph* graph, const std::vector<int>& interestVertices)
		:randomSource(randomSource), graph(graph), interestVertices(interestVertices), n(0), estimateFirstMoment(0), estimateSecondMoment(0), warnedStability(false)
	{}
	void conditionalPMC(ConditionalTurnipInput& input)
	{
		const context::internalGraph& graph = *input.graph;
		const std::vector<int>& interestVertices = input.interestVertices;
		const std::size_t nEdges = boost::num_edges(graph);
		const std::size_t nVertices = boost::num_vertices(graph);
		if (input.edges.size() != nEdges && input.edges.size() != 0)
		{
			throw std::runtime_error("input.edges had an invalid size");
		}
		if (input.edges.size() == 0)
		{
			input.edges.resize(nEdges);
			context::internalGraph::edge_iterator current, end;
			boost::tie(current, end) = boost::edges(graph);
			for (; current != end; current++)
			{
				input.edges[boost::get(boost::edge_index, graph, *current)] = std::make_pair((int)current->m_source, (int)current->m_target);
			}
		}
		if (input.edgeCounts.size() != nEdges)
		{
			throw std::runtime_error("input.edgeCounts had an invalid size");
		}
		//If the vertices are already in the same component, return probability 0. 
		int firstVertex = input.interestVertices[0];
		bool allSameVertex = true;
		for (std::vector<int>::const_iterator i = input.interestVertices.begin(); i != input.interestVertices.end(); i++)
		{
			if (*i != firstVertex)
			{
				allSameVertex = false;
				break;
			}
		}
		if (allSameVertex)
		{
			input.estimateFirstMoment = input.estimateSecondMoment = 0;
			return;
		}
		//If they're never going to be in the same component, return probability 1. 
		if (nEdges == 0)
		{
			input.estimateFirstMoment = input.estimateSecondMoment = 1;
			return;
		}
		boost::random_number_generator<boost::mt19937> generator(input.randomSource);

		input.workingEdgeCounts.clear();
		int originalEdges = 0;
		for (std::vector<int>::iterator i = input.edgeCounts.begin(); i != input.edgeCounts.end(); i++)
		{
			input.workingEdgeCounts.insert(input.workingEdgeCounts.begin(), *i, (int)std::distance(input.edgeCounts.begin(), i));
			originalEdges += *i;
		}
		input.exponentialRates.clear();
		input.exponentialRates.resize(nEdges);
		//The sum of all the conditional probabilities (Used to get estimate and estimate of variance of estimate)
		mpfr_class sumConditional = 0, sumSquaredConditional = 0;
		//This stores the rates that go into the matrix exponential computation
		std::vector<mpfr_class> ratesForPMC, secondRatesForPMC;

		//The initial rate at the start
		mpfr_class sumAllRates = originalEdges * input.exponentialRate;//(originalEdges - input.minimumInoperative) * input.exponentialRate;

		//Only warn about stability once
		input.warnedStability = false;

		//scratch data used to compute conditional probabilities
		std::vector<mpfr_class> computeConditionalProbScratch;
		//Used to hold the connected component IDS
		std::vector<int> componentIDs(nVertices);
		for (std::size_t i = 0; i < input.n; i++)
		{
			//determine which edges are going to be left out (it'll be the last ones of this shuffled vector)
			boost::random_shuffle(input.workingEdgeCounts, generator);
			for (std::size_t j = 0; j < nEdges; j++)
			{
				input.exponentialRates[j] = input.edgeCounts[j] * input.exponentialRate;
			}
			//No edges have yet been seen
			//The first rate is going to be this
			mpfr_class currentRate = sumAllRates;
			//which edge in the permutation are we currently looking at?
			std::size_t permutationCounter = 0;
			//these are going to be the rates for the matrix exponential
			ratesForPMC.clear();
			secondRatesForPMC.clear();
			//reset connected component IDs
			componentIDs.clear();
			componentIDs.insert(componentIDs.begin(), boost::counting_iterator<int>(0), boost::counting_iterator<int>((int)nVertices));
			bool connected = false;
			mpfr_class additionalPart;
			while (true)
			{
				//If we get to a situtation with permutationCounter == input.edgeOrdering.size(), then the graph is 
				//disconnected EVEN when all the edges are repaired. 
				if (permutationCounter == input.workingEdgeCounts.size() - input.minimumInoperative)
				{
					break;
				}
				//get out the edge index
				int edgeIndex = input.workingEdgeCounts[permutationCounter];
				//Get the old component IDs
				int firstComponentID = componentIDs[input.edges[edgeIndex].first];
				int secondComponentID = componentIDs[input.edges[edgeIndex].second];
				//If not, add the current rate
				if (!connected)
				{
					ratesForPMC.push_back(currentRate);
				}
				else secondRatesForPMC.push_back(currentRate);
				//subtract the rate from the current rate
				currentRate -= input.exponentialRate;
				std::replace(componentIDs.begin(), componentIDs.end(), std::max(firstComponentID, secondComponentID), std::min(firstComponentID, secondComponentID));
				//determine whether or not we've hit the critical threshold
				if (!connected)
				{
					connected = true;
					for (std::size_t k = 1; k < interestVertices.size(); k++)
					{
						connected &= (componentIDs[interestVertices[k]] == componentIDs[interestVertices[0]]);
					}
					if (connected)
					{
						additionalPart = computeConditionalProb(ratesForPMC, computeConditionalProbScratch);
						//mpfr_class additionalPart2 = computeConditionalProb(ratesForPMC);
						if ((additionalPart > 1 || additionalPart < 0) && !input.warnedStability)
						{
							input.warnedStability = true;
						}
					}
				}
				permutationCounter++;
			}
			if (!connected)
			{
				additionalPart = 1;
			}
			else if (input.minimumInoperative != 0)
			{
				secondRatesForPMC.push_back(currentRate);
				mpfr_class additionalPart2 = computeConditionalProbLarger(ratesForPMC, secondRatesForPMC, computeConditionalProbScratch);
				additionalPart /= additionalPart2;
			}

			sumConditional += additionalPart;
			sumSquaredConditional += additionalPart*additionalPart;
		}
		input.estimateFirstMoment = sumConditional / input.n;
		input.estimateSecondMoment = sumSquaredConditional / input.n;
		std::string tmp = input.estimateFirstMoment.str();
		tmp = input.estimateSecondMoment.str();
	}
}
