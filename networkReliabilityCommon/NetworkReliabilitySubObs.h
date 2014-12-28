#ifndef NETWORK_RELIABILITY_SUB_OBS_HEADER_GUARD
#define NETWORK_RELIABILITY_SUB_OBS_HEADER_GUARD
#include <boost/random/mersenne_twister.hpp>
#include "Context.h"
#include "EdgeState.h"
#include <boost/shared_array.hpp>
#include "includeMPFR.h"
#include "graphAlgorithms.h"
namespace networkReliability
{
	class NetworkReliabilityObs;
	class NetworkReliabilitySubObs : public boost::noncopyable
	{
	public:
		typedef mpfr_class conditioning_type;
		NetworkReliabilitySubObs(NetworkReliabilitySubObs&& other);
		NetworkReliabilitySubObs(Context const& context, boost::shared_array<EdgeState> state, int radius, int conditioningCount, conditioning_type conditiniongProb);
		NetworkReliabilityObs getObservation(boost::mt19937& randomSource) const;
		const EdgeState* getState() const;
		int getMinCut() const;
		NetworkReliabilitySubObs& operator=(NetworkReliabilitySubObs&& other);
		int getConditioningCount() const;
		const conditioning_type& getConditioningProb() const;
		void getRadius1ReducedGraph(Context::internalGraph& outputGraph, int& minimumInoperative, std::vector<int>& edgeCounts, std::vector<int>& components, boost::detail::depth_first_visit_restricted_impl_helper<Context::internalGraph>::stackType& stack, std::vector<boost::default_color_type>& colorMap) const;
		void getRadius1ReducedGraphNoSelf(Context::internalGraph& outputGraph, int& minimumInoperative, std::vector<int>& edgeCounts, std::vector<int>& components, boost::detail::depth_first_visit_restricted_impl_helper<Context::internalGraph>::stackType& stack, std::vector<boost::default_color_type>& colorMap) const;
		NetworkReliabilitySubObs::conditioning_type getGeneratedObservationConditioningProb() const;
		NetworkReliabilitySubObs copyWithConditioningProb(const conditioning_type& conditioningProb) const;
	private:
		NetworkReliabilitySubObs(Context const& context);
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
