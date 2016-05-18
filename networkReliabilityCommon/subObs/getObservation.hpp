#ifndef NETWORK_RELIABILITY_GET_OBSERVATION_HEADER_GUARD
#define NETWORK_RELIABILITY_GET_OBSERVATION_HEADER_GUARD
#include <boost/noncopyable.hpp>
#include "context.h"
#include "edgeState.h"
namespace networkReliability
{
	namespace subObs
	{
		template<class T> class getObservation : public boost::noncopyable
		{
			public:
				static typename T::observationType get(const T& input, boost::random::mt19937& randomSource)
				{
					const Context& context = input.getContext();
					std::size_t nEdges = context.getNEdges();
					boost::shared_array<edgeState> state(new edgeState[nEdges]);
					typename T::observationConstructorType otherData;
					input.getObservation(state.get(), randomSource, otherData);
					typename T::observationType returnVal(input.getContext(), state, otherData);
					return returnVal;
				}
			private:
				getObservation();
		};
	}
}
#endif
