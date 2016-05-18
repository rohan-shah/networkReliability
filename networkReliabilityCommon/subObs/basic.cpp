#include "subObs/basic.h"
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/math/distributions/binomial.hpp>
#include "graphAlgorithms.h"
#include "seriesParallelReduction.hpp"
namespace networkReliability
{
	namespace subObs
	{
		basic::basic(context const& contextObj, boost::shared_array<edgeState> state, double radius)
			: ::networkReliability::subObs::subObs(contextObj, state, radius), fixedInop(0)
		{
		}
		void basic::initialise()
		{
			potentiallyDisconnected = isSingleComponent(contextObj, state.get(), contextObj.capacityVector, contextObj.depthFirstStack, contextObj.colorVector);
			const std::size_t nEdges = contextObj.getNEdges();
			//Profiling has indicated that the push_back below costs in terms of calls to new[]. Hence the call to reserve first.  
			int couldBeDeactivatedCounter = 0;
			for(std::size_t i = 0; i < nEdges; i++)
			{
				if(state[i] & (UNFIXED_INOP | UNFIXED_OP)) couldBeDeactivatedCounter++;
			}
			couldBeDeactivated.reserve(couldBeDeactivatedCounter);
			for(std::size_t i = 0; i < nEdges; i++)
			{
				if(state[i] & (UNFIXED_INOP | UNFIXED_OP))
				{
					couldBeDeactivated.push_back((int)i);
				}
				else if(state[i] == FIXED_INOP)
				{
					fixedInop++;
				}
			}
		}
		basic::basic(basic&& other)
			: ::networkReliability::subObs::subObs(static_cast<::networkReliability::subObs::subObs&&>(other))
		{
			potentiallyDisconnected = other.potentiallyDisconnected;
			fixedInop = other.fixedInop;
			couldBeDeactivated.swap(other.couldBeDeactivated);
		}
		basic& basic::operator=(basic&& other)
		{
			*static_cast<::networkReliability::subObs::subObs*>(this) = static_cast<::networkReliability::subObs::subObs&&>(other);
			potentiallyDisconnected = other.potentiallyDisconnected;
			fixedInop = other.fixedInop;
			couldBeDeactivated.swap(other.couldBeDeactivated);
			return *this;
		}
		void basic::getObservation(edgeState* newState, boost::mt19937& randomSource, observationConstructorType& otherData) const
		{
			const std::size_t nEdges = contextObj.getNEdges();
			const std::size_t minCut = contextObj.getMinCutEdges();
			memcpy(newState, state.get(), sizeof(edgeState)*nEdges);
			if(radius == 0)
			{
				return;
			}
			const ::TruncatedBinomialDistribution::TruncatedBinomialDistribution& distribution = contextObj.getInopDistribution(std::max(0, (int)minCut - (int)fixedInop), couldBeDeactivated.size(), couldBeDeactivated.size());
			std::size_t nDeactivated = distribution(randomSource);

			for(std::size_t i = 0; i < couldBeDeactivated.size(); i++)
			{
				newState[couldBeDeactivated[i]] = UNFIXED_OP;
			}
			std::vector<int> indices(boost::counting_iterator<int>(0), boost::counting_iterator<int>((int)couldBeDeactivated.size()));
			for(std::size_t i = 0; i < nDeactivated; i++)
			{
				boost::random::uniform_int_distribution<std::size_t> dist(0, indices.size()-1);
				std::size_t generated = dist(randomSource);
				newState[couldBeDeactivated[indices[generated]]] = UNFIXED_INOP;
				std::swap(indices[generated], *indices.rbegin());
				indices.pop_back();
			}
		}
		basic::basic(context const& contextObj, boost::shared_array<edgeState> state, double radius, ::networkReliability::subObs::basicConstructorType&)
			: ::networkReliability::subObs::subObs(contextObj, state, radius), fixedInop(0)
		{
			initialise();
		}
		int basic::getFixedInopCount() const
		{
			return fixedInop;
		}
	}
}
