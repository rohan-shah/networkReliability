#ifndef BINARY_DATASET_HEADER_GUARD
#define BINARY_DATASET_HEADER_GUARD
#include "EdgeState.h"
#include <vector>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/access.hpp>
namespace networkReliability
{
	class binaryDataSet
	{
	protected:
		friend class boost::serialization::access;
		binaryDataSet()
			:nStoredBits(0),storedBits(0)
		{}
		BOOST_SERIALIZATION_SPLIT_MEMBER()
		template<class Archive> void save(Archive & ar, const unsigned int version) const
		{
			ar << data << nStoredBits << storedBits;
		}
		template<class Archive> void load(Archive & ar, const unsigned int version)
		{
			ar >> data >> nStoredBits >> storedBits;
		}
		std::vector<int> data;
		int nStoredBits;
		int storedBits;
	};
	class binaryDataSet1 : public binaryDataSet
	{
	protected:
		binaryDataSet1()
		{}
		friend class boost::serialization::access;
		void add(const EdgeState* state, const std::size_t size);
		BOOST_SERIALIZATION_SPLIT_MEMBER()
		template<class Archive> void save(Archive & ar, const unsigned int version) const
		{
			ar << *static_cast<const binaryDataSet*>(this);
		}
		template<class Archive> void load(Archive & ar, const unsigned int version)
		{
			ar >> *static_cast<binaryDataSet*>(this);
		}
		void expand(int index, std::vector<int>& output);
	};
}
#endif
