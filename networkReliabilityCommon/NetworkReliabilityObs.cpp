#include "NetworkReliabilityObs.h"
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/math/distributions/binomial.hpp>
namespace networkReliability
{
	NetworkReliabilityObs::NetworkReliabilityObs(Context const& context, boost::mt19937& randomSource)
		:context(context), conditioningCount(0), conditioningProb(1)
	{
		const Context::internalGraph& graph = context.getGraph();
		const std::size_t nEdges = boost::num_edges(graph);
		boost::shared_array<EdgeState> state(new EdgeState[nEdges]);

		boost::random::bernoulli_distribution<double> edgeDist(context.getOperationalProbability().toDouble());

		for(std::size_t i = 0; i < nEdges; i++)
		{
			if(edgeDist(randomSource))
			{
				state[i] = UNFIXED_OP;
			}
			else state[i] = UNFIXED_INOP;
		}
		this->state = state;
	}
	NetworkReliabilityObs NetworkReliabilityObs::constructConditional(Context const& context, boost::mt19937& randomSource)
	{
		const Context::internalGraph& graph = context.getGraph();
		const std::size_t nEdges = boost::num_edges(graph);
		const std::size_t minCutEdges = context.getMinCutEdges();
		const ::TruncatedBinomialDistribution::TruncatedBinomialDistribution& dist = context.getDistribution(minCutEdges, nEdges, nEdges);

		const std::size_t nRemovedEdges = dist(randomSource);
		boost::shared_array<EdgeState> state(new EdgeState[nEdges]);
		std::fill(state.get(), state.get() + nEdges, UNFIXED_OP);

		std::vector<int> indices(boost::counting_iterator<int>(0), boost::counting_iterator<int>((int)nEdges));
		for(int i = 0; i < nRemovedEdges; i++)
		{
			boost::random::uniform_int_distribution<int> removedEdgeIndexDistribution(0, (int)indices.size()-1);
			int index = removedEdgeIndexDistribution(randomSource);
			state[indices[index]] = UNFIXED_INOP;
			std::swap(indices[index], *indices.rbegin());
			indices.pop_back();
		}
		boost::math::binomial_distribution<> relevantDistribution((double)nEdges, context.getInoperationalProbabilityD());
		double conditioningProb = 1-boost::math::cdf(relevantDistribution, minCutEdges-1);
		return NetworkReliabilityObs(context, state, (int)minCutEdges, conditioningProb);
	}
	NetworkReliabilityObs::NetworkReliabilityObs(Context const& context, boost::shared_array<EdgeState> state, int conditioningCount, conditioning_type conditioningProb)
		:context(context), state(state), conditioningCount(conditioningCount), conditioningProb(conditioningProb)
	{}
	const EdgeState* NetworkReliabilityObs::getState() const
	{
		return state.get();
	}
	NetworkReliabilityObs& NetworkReliabilityObs::operator=(const NetworkReliabilityObs& other)
	{
		if(&context != &other.context)
		{
			throw std::runtime_error("Internal error");
		}
		state = other.state;
		conditioningCount = other.conditioningCount;
		conditioningProb = other.conditioningProb;
		return *this;
	}
	NetworkReliabilitySubObs NetworkReliabilityObs::getSubObservation(int radius) const
	{
		const Context::internalGraph& graph = context.getGraph();
		std::size_t nEdges = boost::num_edges(graph);
		
		boost::shared_array<EdgeState> newEdgeStates(new EdgeState[nEdges]);
		EdgeState* newEdgeStatesPtr = newEdgeStates.get();
		std::fill(newEdgeStatesPtr, newEdgeStatesPtr + nEdges, FIXED_OP);
		
		const EdgeState* oldEdgeStatesPtr = state.get();

		const int* edgeDistances = context.getEdgeDistances();

		std::size_t sourceEdge = 0;
		while(sourceEdge < nEdges)
		{
			//is this vertex marked as on, for one reason or another? If so continue from here
			if((oldEdgeStatesPtr[sourceEdge] & INOP_MASK) > 0 && newEdgeStatesPtr[sourceEdge] == FIXED_OP)
			{
				newEdgeStatesPtr[sourceEdge] = NEW_FIXED_INOP;

				//Do we find another vertex in our search that is marked on, and is far enough away from the source?
				//If so retain it, it will be our new starting point. 
				//If no such found, we'll continue from finalSearchVertex+1
				bool found = false;
				std::size_t nextSourceEdge = -1;
				//keep copy of source vertex
				std::size_t copiedSourceEdge = sourceEdge;
				//we want to begin on the NEXT vertex
				sourceEdge++;
				while(sourceEdge < nEdges) 
				{
					EdgeState previousState = oldEdgeStatesPtr[sourceEdge];
					if(edgeDistances[copiedSourceEdge + nEdges * sourceEdge] <= radius)
					{
						if(previousState == NEW_FIXED_INOP) newEdgeStatesPtr[sourceEdge] = FIXED_INOP;
						else if(previousState & FIXED_MASK) newEdgeStatesPtr[sourceEdge] = previousState;
						else newEdgeStatesPtr[sourceEdge] = UNFIXED_INOP;
					}
					else if(!found && (previousState & INOP_MASK) > 0 && newEdgeStatesPtr[sourceEdge] == FIXED_OP)
					{
						nextSourceEdge = sourceEdge;
						found = true;
					}
					sourceEdge++;
				}
				//if we found another vertex, continue from there. If no, we're already at finalSearchVertex+1. 
				//Which is where we want to be.
				if(found)
				{
					sourceEdge = nextSourceEdge;
				}
			}
			else sourceEdge++;
		}
		return NetworkReliabilitySubObs(context, newEdgeStates, radius, conditioningCount, conditioningProb);
	}
	NetworkReliabilityObs::NetworkReliabilityObs(NetworkReliabilityObs&& other)
		:context(other.context), state(other.state), conditioningCount(other.conditioningCount), conditioningProb(other.conditioningProb)
	{}
	NetworkReliabilityObs::conditioning_type NetworkReliabilityObs::getConditioningProb() const
	{
		return conditioningProb;
	}
	int NetworkReliabilityObs::getConditioningCount() const
	{
		return conditioningCount;
	}
}