#ifndef NETWORK_REALIABILITY_OBS_HEADER_GUARD
#define NETWORK_REALIABILITY_OBS_HEADER_GUARD
#include "Context.h"
#include "EdgeState.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/shared_array.hpp>
#include "NetworkReliabilitySubObs.h"
namespace networkReliability
{
	class NetworkReliabilityObs
	{
	public:
		typedef mpfr_class conditioning_type;
		NetworkReliabilityObs(Context const& context, boost::mt19937& randomSource);
		NetworkReliabilityObs(Context const& context, boost::shared_array<EdgeState> state, int conditioningCount, conditioning_type conditioningProb);
		NetworkReliabilityObs(NetworkReliabilityObs&& other);
		const EdgeState* getState() const;
		NetworkReliabilityObs& operator=(const NetworkReliabilityObs& other);
		NetworkReliabilitySubObs getSubObservation(double radius) const;
		//construct an object that has some number of edges inoperative - The minimum number needed to disconnect the graph, in fact
		static NetworkReliabilityObs constructConditional(Context const& context, boost::mt19937& randomSource);
		int getConditioningCount() const;
		conditioning_type getConditioningProb() const;
		void getPotentiallyFixed(std::vector<int>& potentiallyFixedIndices, double oldThreshold, double newThreshold, EdgeState* workingMemory) const;
	private:
		Context const& context;
		boost::shared_array<EdgeState> state;
		int conditioningCount;
		conditioning_type conditioningProb;
	};
}
#endif
