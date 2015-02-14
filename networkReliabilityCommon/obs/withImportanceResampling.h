#ifndef NETWORK_RELIABILITY_OBS_WITH_IMPORTANCE_RESAMPLING_HEADER_GUARD
#define NETWORK_RELIABILITY_OBS_WITH_IMPORTANCE_RESAMPLING_HEADER_GUARD
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
		class withImportanceResampling;
	}
	namespace obs
	{
		class withImportanceResampling : public ::networkReliability::withSub
		{
		public:
			template<class T> friend class ::networkReliability::obs::getSubObservation;
			template<class T> friend class ::networkReliability::subObs::getObservation;

			typedef ::networkReliability::subObs::withImportanceResampling subObservationType;
			typedef ::networkReliability::subObs::withImportanceResamplingConstructorType subObservationConstructorType;
			withImportanceResampling(Context const& context, boost::mt19937& randomSource);
			withImportanceResampling(Context const& context, boost::shared_array<EdgeState> state, int conditioningCount, mpfr_class conditioningProb);
			withImportanceResampling(withImportanceResampling&& other);
			withImportanceResampling& operator=(withImportanceResampling&& other);
			//construct an object that has some number of edges inoperative - The minimum number needed to disconnect the graph, in fact
			static ::networkReliability::obs::withImportanceResampling constructConditional(Context const& context, boost::mt19937& randomSource);
			int getConditioningCount() const;
			mpfr_class getConditioningProb() const;
		private:
			withImportanceResampling(Context const& context, boost::shared_array<EdgeState> state, ::networkReliability::obs::withImportanceResamplingConstructorType&);
			void getSubObservation(EdgeState* newState, double radius, subObservationConstructorType& other) const;
			int conditioningCount;
			mpfr_class conditioningProb;
		};
	}
}
#endif
