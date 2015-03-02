#ifndef NETWORK_RELIABILITY_SUB_OBS_WITH_IMPORTANCE_RESAMPLING_HEADER_GUARD
#define NETWORK_RELIABILITY_SUB_OBS_WITH_IMPORTANCE_RESAMPLING_HEADER_GUARD
#include <boost/random/mersenne_twister.hpp>
#include "Context.h"
#include "serializeGMP.hpp"
#include "EdgeState.h"
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
		class withImportanceResampling;
	}
	namespace subObs
	{
		class withImportanceResampling : public ::networkReliability::subObs::subObs
		{
		public:
			template<class T> friend class ::networkReliability::subObs::getObservation;
			template<class T> friend class ::networkReliability::obs::getSubObservation;
			typedef ::networkReliability::obs::withImportanceResampling observationType;
			typedef ::networkReliability::obs::withImportanceResamplingConstructorType observationConstructorType;

			int getFixedInopCount() const;
			const std::vector<int>& getUnknownStateEdges() const;
			withImportanceResampling(withImportanceResampling&& other);
			int getMinCut() const;
			networkReliability::subObs::withImportanceResampling& operator=(withImportanceResampling&& other);
			int getConditioningCount() const;

			const mpfr_class& getConditioningProb() const;
			const mpfr_class& getGeneratedObservationConditioningProb() const;
			networkReliability::subObs::withImportanceResampling copyWithGeneratedObservationConditioningProb(const mpfr_class& conditioningProb) const;
			const mpfr_class& getResamplingProb() const;
		private:
			void initialise();
			withImportanceResampling(Context const& context, boost::shared_array<EdgeState> state, double radius);
			withImportanceResampling(Context const& context, boost::shared_array<EdgeState> state, double radius, ::networkReliability::subObs::withImportanceResamplingConstructorType&);
			void getObservation(EdgeState* state, boost::mt19937& randomSource, observationConstructorType&) const;
			int minCut;
			std::vector<int> unknownState;
			int conditioningCount;
			int fixedInop;
			mpfr_class conditioningProb, generatedObservationConditioningProb;
			int generatedObservationConditioningCount;

			double nextSmallerRadius;
			std::vector<int> boundaryEdges;
			std::vector<int> importanceSamplingEdges;
			mpfr_class resamplingProb;
		};
	}
}
#endif
