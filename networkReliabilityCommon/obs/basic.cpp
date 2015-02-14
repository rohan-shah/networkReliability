#include "basic.h"
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/math/distributions/binomial.hpp>
namespace networkReliability
{
	namespace obs
	{
		basic::basic(Context const& context, boost::shared_array<EdgeState> state, ::networkReliability::obs::basicConstructorType&)
			: ::networkReliability::withSub(context, state)
		{}
		basic::basic(Context const& context, boost::mt19937& randomSource)
			: ::networkReliability::withSub(context, randomSource)
		{}
		basic basic::constructConditional(Context const& context, boost::mt19937& randomSource)
		{
			const std::size_t nEdges = context.getNEdges();
			boost::shared_array<EdgeState> state(new EdgeState[nEdges]);
			NetworkReliabilityObs::constructConditional(context, randomSource, state.get(), false);
			basic retVal(context, state);
			return retVal;
		}
		basic::basic(Context const& context, boost::shared_array<EdgeState> state)
			: ::networkReliability::withSub(context, state)
		{}
		basic& basic::operator=(basic&& other)
		{
			this->::networkReliability::withSub::operator=(static_cast<::networkReliability::withSub&&>(other));
			return *this;
		}
		basic::basic(basic&& other)
			: ::networkReliability::withSub(static_cast<::networkReliability::withSub&&>(other))
		{}
		void basic::getSubObservation(EdgeState* newState, double radius, subObservationConstructorType& other) const
		{
			::networkReliability::withSub::getSubObservation(radius, newState);
		}
	}
}
