#include "NetworkReliabilitySubObsTree.h"
namespace networkReliability
{
	NetworkReliabilitySubObsTree::NetworkReliabilitySubObsTree(boost::archive::binary_iarchive& ar)
	{
		ar >> *this;
	}
	NetworkReliabilitySubObsTree::NetworkReliabilitySubObsTree(boost::archive::text_iarchive& ar)
	{
		ar >> *this;
	}
	NetworkReliabilitySubObsTree::NetworkReliabilitySubObsTree(Context const* externalContext, unsigned int nLevels, unsigned int reservePerLevel, const std::vector<double>& thresholds)
		:thresholds(thresholds)
	{
		parentData.resize(nLevels);
		potentiallyDisconnected.resize(nLevels);
		for(unsigned int i = 0; i < nLevels; i++)
		{
			NetworkReliabilitySubObsCollection currentLevelData(externalContext, thresholds[i]);
			currentLevelData.reserve(reservePerLevel);
			levelData.push_back(std::move(currentLevelData));
			parentData[i].reserve(reservePerLevel);
			potentiallyDisconnected[i].reserve(reservePerLevel);
		}
	}
	std::size_t NetworkReliabilitySubObsTree::getSampleSize(unsigned int level)
	{
		return levelData[level].getSampleSize();
	}
	void NetworkReliabilitySubObsTree::expand(boost::shared_array<EdgeState> state, unsigned int level, unsigned int index) const
	{
		if(level >= levelData.size())
		{
			throw std::runtime_error("Specified level does not exist in call to NetworkReliabilitySubObsTree::expand");
		}
		if(index >= levelData[level].getSampleSize())
		{
			throw std::runtime_error("Specified observation does not exist in specified level, in call to NetworkReliabilitySubObsTree::expand");
		}
		levelData[level].expand(index, state);
	}
}
