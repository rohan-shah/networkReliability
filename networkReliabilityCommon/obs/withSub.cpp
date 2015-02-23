#include "obs/withSub.h"
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/math/distributions/binomial.hpp>
namespace networkReliability
{
	withSub::withSub(Context const& context, boost::mt19937& randomSource)
		: ::networkReliability::NetworkReliabilityObs(context, randomSource)
	{}
	withSub::withSub(Context const& context, boost::shared_array<EdgeState> state)
		: ::networkReliability::NetworkReliabilityObs(context, state)
	{}
	withSub& withSub::operator=(withSub&& other)
	{
		*static_cast<::networkReliability::NetworkReliabilityObs*>(this) = static_cast<::networkReliability::NetworkReliabilityObs&&>(other);
		return *this;
	}
	void withSub::getSubObservation(double radius, EdgeState* newState) const
	{
		getSubObservation(radius, newState, context, state.get());
	}
	void withSub::getSubObservation(double radius, EdgeState* newState, const Context& context, const EdgeState* oldEdgeStatesPtr)
	{
		const Context::internalGraph& graph = context.getGraph();
		std::size_t nEdges = boost::num_edges(graph);
	
		std::fill(newState, newState + nEdges, FIXED_OP);
		const double* edgeDistances = context.getEdgeDistances();

		std::size_t sourceEdge = 0;
		while(sourceEdge < nEdges)
		{
			//is this vertex marked as on, for one reason or another? If so continue from here
			if((oldEdgeStatesPtr[sourceEdge] & INOP_MASK) > 0 && newState[sourceEdge] == FIXED_OP)
			{
				newState[sourceEdge] = FIXED_INOP;

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
						if(previousState & FIXED_MASK) newState[sourceEdge] = previousState;
						else newState[sourceEdge] = UNFIXED_INOP;
					}
					else if(!found && (previousState & INOP_MASK) > 0 && newState[sourceEdge] == FIXED_OP)
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
	}
	withSub::withSub(withSub&& other)
		: ::networkReliability::NetworkReliabilityObs(static_cast<::networkReliability::NetworkReliabilityObs&&>(other))
	{}
}
