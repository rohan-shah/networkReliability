#ifndef NETWORK_RELIABILITY_SUB_OBS_BASIC_HEADER_GUARD
#define NETWORK_RELIABILITY_SUB_OBS_BASIC_HEADER_GUARD
#include <boost/random/mersenne_twister.hpp>
#include "Context.h"
#include "serializeGMP.hpp"
#include "EdgeState.h"
#include <boost/shared_array.hpp>
#include "includeMPFR.h"
#include "graphAlgorithms.h"
#include "subObs/getObservation.hpp"
#include "obs/getSubObservation.hpp"
#include "subObs/subObs.h"
#include "subObsTypes.h"
namespace networkReliability
{
	namespace obs
	{
		class basic;
	}
	namespace subObs
	{
		class basic : public ::networkReliability::subObs::subObs
		{
		public:
			template<class T> friend class ::networkReliability::subObs::getObservation;
			template<class T> friend class ::networkReliability::obs::getSubObservation;

			typedef ::networkReliability::obs::basic observationType;
			typedef ::networkReliability::obs::basicConstructorType observationConstructorType;

			int getFixedInopCount() const;
			const std::vector<int>& getPotentiallyDeactivated() const;
			basic(basic&& other);
			basic(Context const& context, boost::shared_array<EdgeState> state, double radius);
			int getMinCut() const;
			basic& operator=(basic&& other);
		private:
			void initialise();
			basic(Context const& context, boost::shared_array<EdgeState> state, double radius, ::networkReliability::subObs::basicConstructorType&);
			void getObservation(EdgeState* state, boost::mt19937& randomSource, observationConstructorType&) const;
			bool potentiallyDisconnected;
			int fixedInop;
			std::vector<int> couldBeDeactivated;
		};
	}
}

#endif
