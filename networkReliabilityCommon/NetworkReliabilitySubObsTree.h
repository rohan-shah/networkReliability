#ifndef NETWORK_RELIABILITY_SUBOBS_TREE_HEADER_GUARD
#define NETWORK_RELIABILITY_SUBOBS_TREE_HEADER_GUARD
#include "NetworkReliabilitySubObsCollection.h"
#include "NetworkReliabilitySubObs.h"
#include <stdexcept>
#include "binaryDataSet.h"
#include <boost/serialization/access.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
namespace networkReliability
{
	class NetworkReliabilitySubObsTree : public boost::noncopyable
	{
	public:
		friend class boost::serialization::access;
		//The levels go 0, 1, ..., nLevels - 1, with level 0 being the topmost level of the tree
		NetworkReliabilitySubObsTree(Context const* externalContext, unsigned int nLevels, unsigned int reservePerLevel, const std::vector<double>& thresholds);
		NetworkReliabilitySubObsTree(boost::archive::binary_iarchive& ar);
		NetworkReliabilitySubObsTree(boost::archive::text_iarchive& ar);
		void add(const NetworkReliabilitySubObs& subObs, unsigned int level, unsigned int parentIndex, bool potentiallyDisconnected);
		const Context& getContext() const;
		void expand(boost::shared_array<EdgeState> state, unsigned int level, unsigned int index) const;
		std::size_t getSampleSize(unsigned int level);
	private:
		BOOST_SERIALIZATION_SPLIT_MEMBER()
		template<class Archive> void save(Archive& ar, const unsigned int version) const
		{
			std::string typeString = "networkReliabilitySubObsTree";
			ar << typeString;
			if(containedContext)
			{
				ar << *containedContext.get();
			}
			else ar << *externalContext;
			ar << levelData.size();
			for(std::vector<NetworkReliabilitySubObsCollection>::iterator i = levelData.begin(); i != levelData.end(); i++)
			{
				ar << *i;
			}
			ar << parentData;
			ar << potentiallyDisconnected;
			ar << thresholds;
			typeString = "networkReliabilitySubObsTree_end";
			ar << typeString;
		}
		template<class Archive> void load(Archive& ar, const unsigned int version)
		{
			std::string typeString;
			ar >> typeString;
			if(typeString != "networkReliabilitySubObsTree")
			{
				throw std::runtime_error("Incorrect type specifier");
			}
			containedContext.reset(new Context(ar));
			std::size_t levelDataSize;
			ar >> levelDataSize;
			for(std::size_t counter = 0; counter < levelDataSize; counter++)
			{
				NetworkReliabilitySubObsCollection newCollection(ar);
				levelData.push_back(std::move(newCollection));
			}
			ar >> parentData;
			ar >> potentiallyDisconnected;
			ar >> thresholds;
			ar >> typeString;
			if(typeString != "networkReliabilitySubObsTree_end")
			{
				throw std::runtime_error("Incorrect type specifier");
			}
			//Check that vectors all have the same sizes
			std::size_t nLevels = levelData.size();
			if(nLevels != parentData.size() || nLevels != potentiallyDisconnected.size())
			{
				throw std::runtime_error("Inconsistent number of levels");
			}
			for(unsigned int i = 0; i < nLevels; i++)
			{
				std::size_t nObjects = levelData[i].getSampleSize();
				if(parentData[i].size() != nObjects || potentiallyDisconnected[i].size() != nObjects)
				{
					throw std::runtime_error("Inconsistent number of objects at some level");
				}
			}
		}
		//Binary encoding for every observation, at every level of the tree
		std::vector<NetworkReliabilitySubObsCollection> levelData;
		//Values giving the parents of every observation, at every level of the tree
		std::vector<std::vector<int>> parentData;
		std::vector<std::vector<bool>> potentiallyDisconnected;
		std::shared_ptr<Context> containedContext;
		Context const* externalContext;
		std::vector<double> thresholds;
	};
}
#endif
