#ifndef NETWORK_RELIABILITY_SUBOBS_COLLECTION_HEADER_GUARD
#define NETWORK_RELIABILITY_SUBOBS_COLLECTION_HEADER_GUARD
#include "empiricalDistribution.h"
#include "networkReliabilityObs.h"
#include "binaryDataSet.h"
#include <boost/serialization/access.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
namespace networkReliability
{
	class NetworkReliabilityObsCollection : protected binaryDataSet2, public boost::noncopyable
	{
	public:
		friend class boost::serialization::access;
		NetworkReliabilityObsCollection(Context const* externalContext, double radius);
		NetworkReliabilityObsCollection(boost::archive::binary_iarchive& ar);
		NetworkReliabilityObsCollection(boost::archive::text_iarchive& ar);
		NetworkReliabilityObsCollection(NetworkReliabilityObsCollection&& other);
		NetworkReliabilityObsCollection& operator=(NetworkReliabilityObsCollection&& other);
		NetworkReliabilityObsCollection(const empiricalDistribution& other);
		void add(const NetworkReliabilityObs& subObs);
		const Context& getContext() const;
		void expand(int count, boost::shared_array<edgeState> state) const;
		double getRadius() const;
		std::size_t getSampleSize() const;
		void reserve(std::size_t count);
	private:
		BOOST_SERIALIZATION_SPLIT_MEMBER()
		template<class Archive> void save(Archive& ar, const unsigned int version) const
		{
			std::string typeString = "networkReliabilityObsCollection";
			ar << typeString;
			ar << sampleSize;
			ar << radius;
			if(containedContext)
			{
				ar << *containedContext.get();
			}
			else ar << *externalContext;
			ar << *static_cast<const binaryDataSet2*>(this);
			typeString = "networkReliabilityObsCollection_end";
			ar << typeString;
		}
		template<class Archive> void load(Archive& ar, const unsigned int version)
		{
			std::string typeString;
			ar >> typeString;
			if(typeString != "networkReliabilityObsCollection")
			{
				throw std::runtime_error("Incorrect type specifier");
			}
			ar >> sampleSize;
			ar >> radius;
			containedContext.reset(new Context(ar));
			ar >> *static_cast<binaryDataSet2*>(this);
			ar >> typeString;
			if(typeString != "networkReliabilityObsCollection_end")
			{
				throw std::runtime_error("Incorrect type specifier");
			}
		}
		std::size_t sampleSize;
		std::shared_ptr<Context> containedContext;
		Context const* externalContext;
		double radius;
	};
}
#endif
