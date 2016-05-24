#ifndef NETWORK_RELIABILITY_OBS_WITH_RESAMPLING_HEADER_GUARD
#define NETWORK_RELIABILITY_OBS_WITH_RESAMPLING_HEADER_GUARD
#include "context.h"
#include "edgeState.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/shared_array.hpp>
#include "obs/withSub.h"
#include "includeMPFRNetworkReliability.h"
#include "obs/withSub.h"
#include "obs/getSubObservation.hpp"
#include "subObs/getObservation.hpp"
#include "subObsTypes.h"
namespace networkReliability
{
	namespace subObs
	{
		class withResampling;
	}
	namespace obs
	{
		class withResampling : public ::networkReliability::withSub
		{
		public:
			template<class T> friend class ::networkReliability::obs::getSubObservation;
			template<class T> friend class ::networkReliability::subObs::getObservation;

			typedef ::networkReliability::subObs::withResampling subObservationType;
			typedef ::networkReliability::subObs::withResamplingConstructorType subObservationConstructorType;
			withResampling(context const& contextObj, boost::mt19937& randomSource);
			withResampling(context const& contextObj, boost::shared_array<edgeState> state, int conditioningCount, mpfr_class conditioningProb);
			withResampling(withResampling&& other);
			withResampling& operator=(withResampling&& other);
			//construct an object that has some number of edges inoperative - The minimum number needed to disconnect the graph, in fact
			static ::networkReliability::obs::withResampling constructConditional(context const& contextObj, boost::mt19937& randomSource);
			int getConditioningCount() const;
			const mpfr_class& getConditioningProb() const;
		private:
			withResampling(context const& contextObj, boost::shared_array<edgeState> state, ::networkReliability::obs::withResamplingConstructorType&);
			void getSubObservation(edgeState* newState, double radius, subObservationConstructorType& other) const;
			int conditioningCount;
			mpfr_class conditioningProb;
		};
	}
}
#endif
