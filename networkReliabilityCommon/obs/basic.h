#ifndef NETWORK_REALIABILITY_OBS_BASIC_HEADER_GUARD
#define NETWORK_REALIABILITY_OBS_BASIC_HEADER_GUARD
#include "Context.h"
#include "EdgeState.h"
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

			basic(Context const& context, boost::mt19937& randomSource);
			basic(Context const& context, boost::shared_array<EdgeState> state);
			basic(basic&& other);
			basic& operator=(basic&& other);
			static basic constructConditional(Context const& context, boost::mt19937& randomSource);
		private:
			void getSubObservation(EdgeState* newState, double radius, subObservationConstructorType& other) const;
			basic(Context const& context, boost::shared_array<EdgeState> state, ::networkReliability::obs::basicConstructorType&);
		};
	}
}
#endif
