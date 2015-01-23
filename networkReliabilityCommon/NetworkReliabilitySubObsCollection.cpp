#include "NetworkReliabilitySubObsCollection.h"
namespace networkReliability
{
	NetworkReliabilitySubObsCollection::NetworkReliabilitySubObsCollection(Context const* externalContext, double radius)
		:sampleSize(0), externalContext(externalContext), radius(radius)
	{}
	NetworkReliabilitySubObsCollection::NetworkReliabilitySubObsCollection(boost::archive::binary_iarchive& ar)
		:sampleSize(0), externalContext(NULL)
	{
		ar >> *this;	
	}
	NetworkReliabilitySubObsCollection::NetworkReliabilitySubObsCollection(boost::archive::text_iarchive& ar)
		:sampleSize(0), externalContext(NULL)
	{
		ar >> *this;
	}
	void NetworkReliabilitySubObsCollection::reserve(unsigned int count)
	{
		const Context* currentContext;
		if(externalContext) currentContext = externalContext;
		else currentContext = &*containedContext;
		binaryDataSet2::reserve(count * currentContext->getNEdges());
	}
	void NetworkReliabilitySubObsCollection::expand(int count, boost::shared_array<EdgeState> state) const
	{
		const Context* currentContext;
		if(externalContext) currentContext = externalContext;
		else currentContext = &*containedContext;

		binaryDataSet2::expand(count, state.get(), currentContext->getNEdges());
	}
	void NetworkReliabilitySubObsCollection::add(const NetworkReliabilitySubObs& subObs)
	{
		if(subObs.getRadius() != radius)
		{
			throw std::runtime_error("All subObservation objects added to a NetworkReliabilitySubObsCollection must have the same radius");
		}
		//Check that certain key variables are the same - Checking by pointer should be sufficient.
		Context const& subContext = subObs.getContext();
		const Context* currentContext;
		if(externalContext) currentContext = externalContext;
		else currentContext = &*containedContext;

		if(&(subContext.getGraph()) != &(currentContext->getGraph()) || &(subContext.getDirectedGraph()) != &(currentContext->getDirectedGraph()) || subContext.getEdgeDistances() != currentContext->getEdgeDistances() || &(subContext.getOperationalProbability()) != &(currentContext->getOperationalProbability()) || subContext.getNEdges() != currentContext->getNEdges())
		{
			throw std::runtime_error("subObservation object added to NetworkReliabilitySubObsCollection had wrong Context object");
		}
		static_cast<binaryDataSet2*>(this)->add(subObs.getState(), subObs.getContext().getNEdges());
		sampleSize++;
	}
	double NetworkReliabilitySubObsCollection::getRadius() const
	{
		return radius;
	}
	const Context& NetworkReliabilitySubObsCollection::getContext() const
	{
		if(containedContext) return *containedContext.get();
		if(externalContext) return *externalContext;
		throw std::runtime_error("Invalid state for NetworkReliabilitySubObsCollection");
	}
	NetworkReliabilitySubObsCollection::NetworkReliabilitySubObsCollection(NetworkReliabilitySubObsCollection&& other)
		:binaryDataSet2(std::move(other)), sampleSize(other.sampleSize), containedContext(other.containedContext), externalContext(other.externalContext), radius(other.radius)
	{}
	NetworkReliabilitySubObsCollection& NetworkReliabilitySubObsCollection::operator=(NetworkReliabilitySubObsCollection&& other)
	{
		sampleSize = other.sampleSize;
		containedContext = other.containedContext;
		externalContext = other.externalContext;
		radius = other.radius;
		*static_cast<binaryDataSet2*>(this) = std::move(static_cast<binaryDataSet2&&>(other));
		return *this;
	}
	std::size_t NetworkReliabilitySubObsCollection::getSampleSize() const
	{
		return sampleSize;
	}
}
