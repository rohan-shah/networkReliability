#ifndef OBS_GET_SUB_OBSERVATION_HEADER_GUARD
#define OBS_GET_SUB_OBSERVATION_HEADER_GUARD
#include "EdgeState.h"
#include "Context.h"
#include <boost/noncopyable.hpp>
namespace networkReliability
{
	namespace obs
	{
		template<class T> class getSubObservation : public boost::noncopyable
		{
		public:
			static typename T::subObservationType get(const T& input, double radius)
			{
				const Context& context = input.getContext();
				typename T::subObservationConstructorType otherData;
				boost::shared_array<EdgeState> newState(new EdgeState[context.getNEdges()]);
				input.getSubObservation(newState.get(), radius, otherData);

				typename T::subObservationType retVal(context, newState, radius, otherData);
				return retVal;
			}
			template<typename U> static typename T::subObservationType get(const T& input, double radius, const U& aux)
			{
				const Context& context = input.getContext();
				typename T::subObservationConstructorType otherData;
				boost::shared_array<EdgeState> newState(new EdgeState[context.getNEdges()]);
				input.getSubObservation(newState.get(), radius, otherData, aux);

				typename T::subObservationType retVal(context, newState, radius, otherData);
				return retVal;
			}
		private:
			getSubObservation();
		};
	}
}
#endif
