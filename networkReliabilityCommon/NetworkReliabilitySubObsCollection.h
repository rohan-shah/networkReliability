#ifndef NETWORK_RELIABILITY_SUBOBS_COLLECTION_HEADER_GUARD
#define NETWORK_RELIABILITY_SUBOBS_COLLECTION_HEADER_GUARD
#include "empiricalDistribution.h"
#include "NetworkReliabilitySubObs.h"
#include "binaryDataSet.h"
#include <boost/serialization/access.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
namespace networkReliability
{
	class NetworkReliabilitySubObs;
	class NetworkReliabilitySubObsCollection : protected binaryDataSet2, public boost::noncopyable
	{
	public:
		friend class boost::serialization::access;
		NetworkReliabilitySubObsCollection(Context const* externalContext, double radius);
		NetworkReliabilitySubObsCollection(boost::archive::binary_iarchive& ar);
		NetworkReliabilitySubObsCollection(boost::archive::text_iarchive& ar);
		NetworkReliabilitySubObsCollection(NetworkReliabilitySubObsCollection&& other);
		NetworkReliabilitySubObsCollection& operator=(NetworkReliabilitySubObsCollection&& other);
		NetworkReliabilitySubObsCollection(const empiricalDistribution& other);
		void add(const NetworkReliabilitySubObs& subObs);
		const Context& getContext() const;
		void expand(int count, boost::shared_array<EdgeState> state) const;
		double getRadius() const;
		std::size_t getSampleSize() const;
		void reserve(unsigned int count);
	private:
		BOOST_SERIALIZATION_SPLIT_MEMBER()
		template<class Archive> void save(Archive& ar, const unsigned int version) const
		{
			std::string typeString = "networkReliabilitySubObsCollection";
			ar << typeString;
			ar << sampleSize;
			ar << radius;
			if(containedContext)
			{
				ar << *containedContext.get();
			}
			else ar << *externalContext;
			ar << *static_cast<const binaryDataSet2*>(this);
			typeString = "networkReliabilitySubObsCollection_end";
			ar << typeString;
		}
		template<class Archive> void load(Archive& ar, const unsigned int version)
		{
			std::string typeString;
			ar >> typeString;
			if(typeString != "networkReliabilitySubObsCollection")
			{
				throw std::runtime_error("Incorrect type specifier");
			}
			ar >> sampleSize;
			ar >> radius;
			containedContext.reset(new Context(ar));
			ar >> *static_cast<binaryDataSet2*>(this);
			ar >> typeString;
			if(typeString != "networkReliabilitySubObsCollection_end")
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
