#include "obs/withImportanceResampling.h"
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/math/distributions/binomial.hpp>
namespace networkReliability
{
	namespace obs
	{
		withImportanceResampling::withImportanceResampling(Context const& context, boost::mt19937& randomSource)
			: ::networkReliability::withSub(context, randomSource), conditioningCount(0), conditioningProb(1) 
		{}
		withImportanceResampling withImportanceResampling::constructConditional(Context const& context, boost::mt19937& randomSource)
		{
			const Context::internalGraph& graph = context.getGraph();
			const std::size_t nEdges = boost::num_edges(graph);
			const std::size_t minCutEdges = context.getMinCutEdges();
			boost::shared_array<EdgeState> state(new EdgeState[nEdges]);
			::networkReliability::NetworkReliabilityObs::constructConditional(context, randomSource, state.get(), false);

			const ::TruncatedBinomialDistribution::TruncatedBinomialDistribution& originalDist = context.getOpDistribution(0, nEdges, nEdges);
			mpfr_class conditioningProb = originalDist.getCumulativeProbability((int)(nEdges - minCutEdges));
			return withImportanceResampling(context, state, (int)minCutEdges, conditioningProb);
		}
		withImportanceResampling::withImportanceResampling(Context const& context, boost::shared_array<EdgeState> state, int conditioningCount, mpfr_class conditioningProb)
			:  ::networkReliability::withSub(context, state), conditioningCount(conditioningCount), conditioningProb(conditioningProb)
		{}
		withImportanceResampling& withImportanceResampling::operator=(withImportanceResampling&& other)
		{
			this->::networkReliability::withSub::operator=(static_cast<::networkReliability::withSub&&>(other));
			conditioningCount = other.conditioningCount;
			conditioningProb = other.conditioningProb;
			return *this;
		}
		withImportanceResampling::withImportanceResampling(Context const& context, boost::shared_array<EdgeState> state, ::networkReliability::obs::withImportanceResamplingConstructorType& other)
			: ::networkReliability::withSub(context, state), conditioningCount(other.conditioningCount), conditioningProb(other.conditioningProb)
		{
		}
		void withImportanceResampling::getSubObservation(EdgeState* newState, double radius, subObservationConstructorType& other, double nextSmallerRadius) const
		{
			other.nextRadius = nextSmallerRadius;
			other.conditioningCount = conditioningCount;
			other.conditioningProb = conditioningProb;

			const Context::internalGraph& graph = context.getGraph();
			std::size_t nEdges = boost::num_edges(graph);
			
			std::fill(newState, newState + nEdges, FIXED_OP);
			
			const EdgeState* oldEdgeStatesPtr = state.get();
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
						double distance = edgeDistances[copiedSourceEdge + nEdges * sourceEdge];
						if(nextSmallerRadius < distance && distance <= radius)
						{
							if(previousState & FIXED_MASK) newState[sourceEdge] = previousState;
							else if(newState[sourceEdge] != UNFIXED_INOP) newState[sourceEdge] = UNFIXED_OP;
						}
						else if(edgeDistances[copiedSourceEdge + nEdges * sourceEdge] <= nextSmallerRadius)
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
			other.boundaryEdges.clear();
			other.interiorEdges.clear();
			//go through and work out which ones are potentially fixed in the next step
			for(EdgeState* newEdgeStateIterator = newState; newEdgeStateIterator != newState + nEdges; newEdgeStateIterator++)
			{
				if(*newEdgeStateIterator == UNFIXED_OP) other.boundaryEdges.push_back((int)(newEdgeStateIterator - newState));
				else if(*newEdgeStateIterator == UNFIXED_INOP) other.interiorEdges.push_back((int)(newEdgeStateIterator - newState));
			}
		}
		withImportanceResampling::withImportanceResampling(withImportanceResampling&& other)
			: ::networkReliability::withSub(static_cast<withImportanceResampling&&>(other)), conditioningCount(other.conditioningCount), conditioningProb(other.conditioningProb)
		{}
		mpfr_class withImportanceResampling::getConditioningProb() const
		{
			return conditioningProb;
		}
		int withImportanceResampling::getConditioningCount() const
		{
			return conditioningCount;
		}
	}
}
