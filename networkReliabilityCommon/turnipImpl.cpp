#include "turnipImpl.h"
#include <boost/iterator/counting_iterator.hpp>
#include <boost/random/random_number_generator.hpp>
#include <boost/range/algorithm/random_shuffle.hpp>
#include "computeConditionalProb.h"
namespace networkReliability
{
	TurnipEqualRateInput::TurnipEqualRateInput(boost::mt19937& randomSource, const context::internalGraph* graph, const std::vector<int>& interestVertices)
		:randomSource(randomSource), graph(graph), interestVertices(interestVertices), n(0), estimateFirstMoment(0), estimateSecondMoment(0), warnedStability(false)
	{}
	void turnipEqualRate(TurnipEqualRateInput& input)
	{
		const context::internalGraph& graph = *input.graph;
		const std::vector<int>& interestVertices = input.interestVertices;
		const std::size_t nEdges = boost::num_edges(graph);
		const std::size_t nVertices = boost::num_vertices(graph);
		if (input.edges.size() != nEdges && input.edges.size() != 0)
		{
			throw std::runtime_error("input.nEdges must had an invalid size");
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
		//The sum of all the conditional probabilities (Used to get estimate and estimate of variance of estimate)
		mpfr_class sumConditional = 0, sumSquaredConditional = 0;
		//This stores the rates that go into the matrix exponential computation
		std::vector<mpfr_class> ratesForPMC;
		boost::random_number_generator<boost::mt19937> generator(input.randomSource);
		std::vector<int> edgeOrdering(boost::counting_iterator<int>(0), boost::counting_iterator<int>((int)nEdges));

		//The initial rate at the start
		mpfr_class sumAllRates = nEdges*input.exponentialRate;

		//The edges which are still under consideration
		std::vector<bool> alreadySeen(nEdges);
		//Only warn about stability once
		input.warnedStability = false;

		//scratch data used to compute conditional probabilities
		std::vector<mpfr_class> computeConditionalProbScratch;
		//Used to hold the connected component IDS
		std::vector<int> componentIDs(nVertices);
		for (std::size_t i = 0; i < input.n; i++)
		{
			//Simulate permutation
			boost::random_shuffle(edgeOrdering, generator);
			//No edges have yet been seen
			std::fill(alreadySeen.begin(), alreadySeen.end(), false);
			//The first rate is going to be this
			mpfr_class currentRate = sumAllRates;
			//which edge in the permutation are we currently looking at?
			std::size_t permutationCounter = 0;
			//these are going to be the rates for the matrix exponential
			ratesForPMC.clear();
			//reset connected component IDs
			componentIDs.clear();
			componentIDs.insert(componentIDs.begin(), boost::counting_iterator<int>(0), boost::counting_iterator<int>((int)nVertices));
			while (true)
			{
				//If we get to a situtation with permutationCounter == nEdges, then the graph is 
				//disconnected EVEN when all the edges are repaired. 
				if (permutationCounter == nEdges)
				{
					break;
				}
				//get out the edge index
				int edgeIndex = edgeOrdering[permutationCounter];
				//Has this edge been removed from consideration already?
				if (!alreadySeen[edgeIndex])
				{
					//If not, add the current rate
					ratesForPMC.push_back(currentRate);
					//Mark the edge as seen
					alreadySeen[edgeIndex] = true;
					//subtract the rate from the current rate
					currentRate -= input.exponentialRate;
					//Get the old component IDs
					int firstComponentID = componentIDs[input.edges[edgeIndex].first];
					int secondComponentID = componentIDs[input.edges[edgeIndex].second];
					//go through all the edges that haven't been seen yet. 
					for (std::size_t j = 0; j < nEdges; j++)
					{
						if (!alreadySeen[j])
						{
							int firstVertex = input.edges[j].first, secondVertex = input.edges[j].second;
							if ((componentIDs[firstVertex] == secondComponentID || componentIDs[firstVertex] == firstComponentID) && (componentIDs[secondVertex] == firstComponentID || componentIDs[secondVertex] == secondComponentID))
							{
								alreadySeen[j] = true;
								currentRate -= input.exponentialRate;
							}
						}
					}
					std::replace(componentIDs.begin(), componentIDs.end(), std::max(firstComponentID, secondComponentID), std::min(firstComponentID, secondComponentID));
					//determine whether or not we've hit the critical threshold
					bool connected = true;
					for (std::size_t k = 1; k < interestVertices.size(); k++)
					{
						connected &= (componentIDs[interestVertices[k]] == componentIDs[interestVertices[0]]);
					}
					if (connected) break;
				}
				permutationCounter++;
			}
			mpfr_class additionalPart;
			if (permutationCounter == nEdges)
			{
				additionalPart = 1;
			}
			else
			{
				additionalPart = computeConditionalProb(ratesForPMC, computeConditionalProbScratch);
				//mpfr_class additionalPart2 = computeConditionalProb(ratesForPMC);
				if ((additionalPart > 1 || additionalPart < 0) && !input.warnedStability)
				{
					input.warnedStability = true;
				}
			}
			sumConditional += additionalPart;
			sumSquaredConditional += additionalPart*additionalPart;
		}
		input.estimateFirstMoment = sumConditional / input.n;
		input.estimateSecondMoment = sumSquaredConditional / input.n;
	}
}
