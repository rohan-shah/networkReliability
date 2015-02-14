#ifndef NETWORK_RELIABILITY_OBS_WITH_RESAMPLING_HEADER_GUARD
#define NETWORK_RELIABILITY_OBS_WITH_RESAMPLING_HEADER_GUARD
#include "Context.h"
#include "EdgeState.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/shared_array.hpp>
#include "obs/withSub.h"
#include "includeMPFR.h"
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
			withResampling(Context const& context, boost::mt19937& randomSource);
			withResampling(Context const& context, boost::shared_array<EdgeState> state, int conditioningCount, mpfr_class conditioningProb);
			withResampling(withResampling&& other);
			withResampling& operator=(withResampling&& other);
			//construct an object that has some number of edges inoperative - The minimum number needed to disconnect the graph, in fact
			static ::networkReliability::obs::withResampling constructConditional(Context const& context, boost::mt19937& randomSource);
			int getConditioningCount() const;
			const mpfr_class& getConditioningProb() const;
		private:
			withResampling(Context const& context, boost::shared_array<EdgeState> state, ::networkReliability::obs::withResamplingConstructorType&);
			void getSubObservation(EdgeState* newState, double radius, subObservationConstructorType& other) const;
			int conditioningCount;
			mpfr_class conditioningProb;
		};
	}
}
#endif
