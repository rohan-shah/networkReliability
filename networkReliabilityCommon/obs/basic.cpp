#include "basic.h"
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/math/distributions/binomial.hpp>
namespace networkReliability
{
	namespace obs
	{
		basic::basic(context const& contextObj, boost::shared_array<edgeState> state, ::networkReliability::obs::basicConstructorType&)
			: ::networkReliability::withSub(contextObj, state)
		{}
		basic::basic(context const& contextObj, boost::mt19937& randomSource)
			: ::networkReliability::withSub(contextObj, randomSource)
		{}
		basic basic::constructConditional(context const& contextObj, boost::mt19937& randomSource)
		{
			const std::size_t nEdges = contextObj.getNEdges();
			boost::shared_array<edgeState> state(new edgeState[nEdges]);
			NetworkReliabilityObs::constructConditional(contextObj, randomSource, state.get(), false);
			basic retVal(contextObj, state);
			return retVal;
		}
		basic::basic(context const& contextObj, boost::shared_array<edgeState> state)
			: ::networkReliability::withSub(contextObj, state)
		{}
		basic& basic::operator=(basic&& other)
		{
			this->::networkReliability::withSub::operator=(static_cast<::networkReliability::withSub&&>(other));
			return *this;
		}
		basic::basic(basic&& other)
			: ::networkReliability::withSub(static_cast<::networkReliability::withSub&&>(other))
		{}
		void basic::getSubObservation(edgeState* newState, double radius, subObservationConstructorType& other) const
		{
			::networkReliability::withSub::getSubObservation(radius, newState);
		}
	}
}
