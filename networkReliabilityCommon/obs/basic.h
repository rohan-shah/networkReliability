#ifndef NETWORK_REALIABILITY_OBS_BASIC_HEADER_GUARD
#define NETWORK_REALIABILITY_OBS_BASIC_HEADER_GUARD
#include "context.h"
#include "edgeState.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/shared_array.hpp>
#include "obs/withSub.h"
#include "obs/getSubObservation.hpp"
#include "subObs/getObservation.hpp"
#include "subObsTypes.h"
namespace networkReliability
{
	namespace subObs
	{
		class basic;
	}
	namespace obs
	{
		class basic : public ::networkReliability::withSub
		{
		public:
			template<class T> friend class ::networkReliability::obs::getSubObservation;
			template<class T> friend class ::networkReliability::subObs::getObservation;

			typedef ::networkReliability::subObs::basic subObservationType;
			typedef ::networkReliability::subObs::basicConstructorType subObservationConstructorType;

			basic(context const& contextObj, boost::mt19937& randomSource);
			basic(context const& contextObj, boost::shared_array<edgeState> state);
			basic(basic&& other);
			basic& operator=(basic&& other);
			static basic constructConditional(context const& contextObj, boost::mt19937& randomSource);
		private:
			void getSubObservation(edgeState* newState, double radius, subObservationConstructorType& other) const;
			basic(context const& contextObj, boost::shared_array<edgeState> state, ::networkReliability::obs::basicConstructorType&);
		};
	}
}
#endif
