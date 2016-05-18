#ifndef NETWORK_RELIABILITY_OBS_WITHSUB_HEADER_GUARD
#define NETWORK_RELIABILITY_OBS_WITHSUB_HEADER_GUARD
#include "context.h"
#include "edgeState.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/shared_array.hpp>
#include "depth_first_search_restricted.hpp"
#include "networkReliabilityObs.h"
namespace networkReliability
{
	class withSub : public ::networkReliability::NetworkReliabilityObs
	{
	public:
		static void getSubObservation(double radius, edgeState* newState, const Context& context, const edgeState* oldEdgeStatesPtr);
	protected:
		withSub(Context const& context, boost::mt19937& randomSource);
		withSub(Context const& context, boost::shared_array<edgeState> state);
		withSub(withSub&& other);
		withSub& operator=(withSub&& other);
		void getSubObservation(double radius, edgeState* outputState) const;
	};
}
#endif
