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
		std::vector<unsigned int> data;
		int nStoredBits;
		unsigned int storedBits;
	};
	class binaryDataSet1 : public binaryDataSet
	{
	public:
		void add(const EdgeState* state, const std::size_t size);
	protected:
		binaryDataSet1()
		{}
		friend class boost::serialization::access;
		BOOST_SERIALIZATION_SPLIT_MEMBER()
		template<class Archive> void save(Archive & ar, const unsigned int version) const
		{
			ar << *static_cast<const binaryDataSet*>(this);
		}
		template<class Archive> void load(Archive & ar, const unsigned int version)
		{
			ar >> *static_cast<binaryDataSet*>(this);
		}
		void expand(int index, std::vector<int>& output) const;
	};
	class binaryDataSet2 : public binaryDataSet
	{
	public:
		void add(const EdgeState* state, const std::size_t size);
	protected:
		binaryDataSet2()
		{}
		friend class boost::serialization::access;
		BOOST_SERIALIZATION_SPLIT_MEMBER()
		template<class Archive> void save(Archive & ar, const unsigned int version) const
		{
			ar << *static_cast<const binaryDataSet*>(this);
		}
		template<class Archive> void load(Archive & ar, const unsigned int version)
		{
			ar >> *static_cast<binaryDataSet*>(this);
		}
		void expand(int index, EdgeState* output, const std::size_t nEdges) const;
	};
}
#endif
