#ifndef NETWORK_RELIABILITY_SUB_OBS_WITH_RESAMPLING_HEADER_GUARD
#define NETWORK_RELIABILITY_SUB_OBS_WITH_RESAMPLING_HEADER_GUARD
#include <boost/random/mersenne_twister.hpp>
#include "context.h"
#include "serializeGMP.hpp"
#include "edgeState.h"
#include <boost/shared_array.hpp>
#include "includeMPFR.h"
#include "graphAlgorithms.h"
#include "subObs/subObs.h"
#include "subObs/getObservation.hpp"
#include "obs/getSubObservation.hpp"
#include "subObsTypes.h"
namespace networkReliability
{
	namespace obs
	{
		class withResampling;
	}
	namespace subObs
	{
		class withResampling : public ::networkReliability::subObs::subObs
		{
		public:
			template<class T> friend class ::networkReliability::subObs::getObservation;
			template<class T> friend class ::networkReliability::obs::getSubObservation;
			typedef ::networkReliability::obs::withResampling observationType;
			typedef ::networkReliability::obs::withResamplingConstructorType observationConstructorType;

			int getFixedInopCount() const;
			const std::vector<int>& getPotentiallyDeactivated() const;
			withResampling(withResampling&& other);
			withResampling(Context const& context, boost::shared_array<edgeState> state, double radius, int conditioningCount, mpfr_class conditiniongProb);
			int getMinCut() const;
			networkReliability::subObs::withResampling& operator=(withResampling&& other);
			int getConditioningCount() const;

			const mpfr_class& getConditioningProb() const;
			const mpfr_class& getGeneratedObservationConditioningProb() const;
			networkReliability::subObs::withResampling copyWithGeneratedObservationConditioningProb(const mpfr_class& conditioningProb) const;
		private:
			withResampling(Context const& context, boost::shared_array<edgeState> state, double radius, ::networkReliability::subObs::withResamplingConstructorType&);
			withResampling(Context const& context, boost::shared_array<edgeState> state, double radius);
			void initialise();
			void getObservation(edgeState* state, boost::mt19937& randomSource, observationConstructorType&) const;
			int minCut;
			std::vector<int> couldBeDeactivated;
			int conditioningCount;
			int fixedInop;
			mpfr_class conditioningProb, generatedObservationConditioningProb;
			int generatedObservationConditioningCount;
		};
	}
}
#endif
