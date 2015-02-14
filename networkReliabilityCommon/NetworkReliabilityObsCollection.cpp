#include "NetworkReliabilityObsCollection.h"
namespace networkReliability
{
	NetworkReliabilityObsCollection::NetworkReliabilityObsCollection(Context const* externalContext, double radius)
		:sampleSize(0), externalContext(externalContext), radius(radius)
	{}
	NetworkReliabilityObsCollection::NetworkReliabilityObsCollection(boost::archive::binary_iarchive& ar)
		:sampleSize(0), externalContext(NULL)
	{
		ar >> *this;	
	}
	NetworkReliabilityObsCollection::NetworkReliabilityObsCollection(boost::archive::text_iarchive& ar)
		:sampleSize(0), externalContext(NULL)
	{
		ar >> *this;
	}
	NetworkReliabilityObsCollection::NetworkReliabilityObsCollection(const empiricalDistribution& other)
		:sampleSize(other.getNSamples()), radius(0)
	{
		externalContext = &other.getContext();
		std::size_t nEdges = externalContext->getNEdges();
		std::vector<int> vectorStates(nEdges);
		boost::scoped_array<EdgeState> state(new EdgeState[nEdges]);
		for(std::size_t i = 0; i < sampleSize; i++)
		{
			other.expand(i, vectorStates);
			for(std::size_t j = 0; j < nEdges; j++)
			{
				if(vectorStates[j]) state[j] = FIXED_OP;
				else state[j] = FIXED_INOP;
			}
			static_cast<binaryDataSet2*>(this)->add(state.get(), nEdges);
		}
	}
	void NetworkReliabilityObsCollection::reserve(unsigned int count)
	{
		const Context* currentContext;
		if(externalContext) currentContext = externalContext;
		else currentContext = &*containedContext;
		binaryDataSet2::reserve(count * currentContext->getNEdges());
	}
	void NetworkReliabilityObsCollection::expand(int count, boost::shared_array<EdgeState> state) const
	{
		const Context* currentContext;
		if(externalContext) currentContext = externalContext;
		else currentContext = &*containedContext;

		binaryDataSet2::expand(count, state.get(), currentContext->getNEdges());
	}
	void NetworkReliabilityObsCollection::add(const NetworkReliabilityObs& obs)
	{
		//Check that certain key variables are the same - Checking by pointer should be sufficient.
		Context const& obsContext = obs.getContext();
		const Context* currentContext;
		if(externalContext) currentContext = externalContext;
		else currentContext = &*containedContext;

		if(&(obsContext.getGraph()) != &(currentContext->getGraph()) || &(obsContext.getDirectedGraph()) != &(currentContext->getDirectedGraph()) || obsContext.getEdgeDistances() != currentContext->getEdgeDistances() || &(obsContext.getOperationalProbability()) != &(currentContext->getOperationalProbability()) || obsContext.getNEdges() != currentContext->getNEdges())
		{
			throw std::runtime_error("observation object added to NetworkReliabilityObsCollection had wrong Context object");
		}
		static_cast<binaryDataSet2*>(this)->add(obs.getState(), obs.getContext().getNEdges());
		sampleSize++;
	}
	double NetworkReliabilityObsCollection::getRadius() const
	{
		return radius;
	}
	const Context& NetworkReliabilityObsCollection::getContext() const
	{
		if(containedContext) return *containedContext.get();
		if(externalContext) return *externalContext;
		throw std::runtime_error("Invalid state for NetworkReliabilityObsCollection");
	}
	NetworkReliabilityObsCollection::NetworkReliabilityObsCollection(NetworkReliabilityObsCollection&& other)
		:binaryDataSet2(std::move(other)), sampleSize(other.sampleSize), containedContext(other.containedContext), externalContext(other.externalContext), radius(other.radius)
	{}
	NetworkReliabilityObsCollection& NetworkReliabilityObsCollection::operator=(NetworkReliabilityObsCollection&& other)
	{
		sampleSize = other.sampleSize;
		containedContext = other.containedContext;
		externalContext = other.externalContext;
		radius = other.radius;
		*static_cast<binaryDataSet2*>(this) = std::move(static_cast<binaryDataSet2&&>(other));
		return *this;
	}
	std::size_t NetworkReliabilityObsCollection::getSampleSize() const
	{
		return sampleSize;
	}
}
