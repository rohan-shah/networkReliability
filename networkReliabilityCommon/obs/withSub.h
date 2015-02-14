#ifndef NETWORK_RELIABILITY_OBS_WITHSUB_HEADER_GUARD
#define NETWORK_RELIABILITY_OBS_WITHSUB_HEADER_GUARD
#include "Context.h"
#include "EdgeState.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/shared_array.hpp>
#include "depth_first_search_restricted.hpp"
#include "NetworkReliabilityObs.h"
namespace networkReliability
{
	class withSub : public ::networkReliability::NetworkReliabilityObs
	{
	public:
	protected:
		withSub(Context const& context, boost::mt19937& randomSource);
		withSub(Context const& context, boost::shared_array<EdgeState> state);
		withSub(withSub&& other);
		withSub& operator=(withSub&& other);
		void getSubObservation(double radius, EdgeState* outputState) const;
	};
}
#endif
