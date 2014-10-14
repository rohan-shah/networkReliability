#ifndef NETWORK_RELIABILITY_SUB_OBS_HEADER_GUARD
#define NETWORK_RELIABILITY_SUB_OBS_HEADER_GUARD
#include <boost/random/mersenne_twister.hpp>
#include "Context.h"
#include "EdgeState.h"
#include <boost/shared_array.hpp>
#include "includeMPIRXX.h"
namespace networkReliability
{
	class NetworkReliabilityObs;
	class NetworkReliabilitySubObs : public boost::noncopyable
	{
	public:
		typedef mpf_class conditioning_type;
		NetworkReliabilitySubObs(NetworkReliabilitySubObs&& other);
		NetworkReliabilitySubObs(Context const& context, boost::shared_array<EdgeState> state, int radius, int conditioningCount, conditioning_type conditiniongProb);
		NetworkReliabilityObs getObservation(boost::mt19937& randomSource) const;
		const EdgeState* getState() const;
		int getMinCut() const;
		NetworkReliabilitySubObs& operator=(NetworkReliabilitySubObs&& other);
		int getConditioningCount() const;
	private:
		Context const& context;
		boost::shared_array<EdgeState> state;
		int radius;
		int minCut;
		std::vector<int> couldBeDeactivated;
		int conditioningCount;
		int fixedInop;
		conditioning_type conditioningProb;
	};
}
#endif