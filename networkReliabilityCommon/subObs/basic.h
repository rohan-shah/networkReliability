#ifndef NETWORK_RELIABILITY_SUB_OBS_BASIC_HEADER_GUARD
#define NETWORK_RELIABILITY_SUB_OBS_BASIC_HEADER_GUARD
#include <boost/random/mersenne_twister.hpp>
#include "context.h"
#include "serializeGMPNetworkReliability.hpp"
#include "edgeState.h"
#include <boost/shared_array.hpp>
#include "includeMPFRNetworkReliability.h"
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
			basic(context const& contextObj, boost::shared_array<edgeState> state, double radius);
			int getMinCut() const;
			basic& operator=(basic&& other);
		private:
			void initialise();
			basic(context const& contextObj, boost::shared_array<edgeState> state, double radius, ::networkReliability::subObs::basicConstructorType&);
			void getObservation(edgeState* state, boost::mt19937& randomSource, observationConstructorType&) const;
			bool potentiallyDisconnected;
			int fixedInop;
			std::vector<int> couldBeDeactivated;
		};
	}
}

#endif
