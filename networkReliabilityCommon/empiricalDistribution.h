#ifndef EMPIRICAL_DISTRIBUTION_HEADER_GUARD
#define EMPIRICAL_DISTRIBUTION_HEADER_GUARD
#include <boost/noncopyable.hpp>
#include "astCode.hpp"
#include "EdgeState.h"
#include <string>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/vector.hpp>
#include <stdexcept>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
namespace networkReliability
{
	class empiricalDistribution : public boost::noncopyable
	{
	public:
		empiricalDistribution(boost::archive::text_iarchive& ar)
		{
			ar >> *this;
		}
		empiricalDistribution(boost::archive::binary_iarchive& ar)
		{
			ar >> *this;
		}
		friend class boost::serialization::access;
		empiricalDistribution(empiricalDistribution&& other);
		empiricalDistribution& operator=(empiricalDistribution&& other);
		empiricalDistribution(bool isWeighted, std::size_t nEdges);
		void hintDataCount(std::size_t size);
		bool save(std::string file);
		void add(const EdgeState* state);
		void add(const EdgeState* state, double weight);
		void expand(int count, std::vector<int>& output);
		std::size_t getNSamples() const;
		std::size_t getNEdges() const;
		bool isWeighted() const;
		double getWeight(std::size_t index) const;
	private:
		template<class Archive> void save(Archive & ar, const unsigned int version) const
		{
			std::string typeString = "empiricalDistribution";
			ar << typeString;
			if(_isWeighted)
			{
				std::string weighted = "weighted";
				ar << weighted;
			}
			else
			{
				std::string unweighted = "unweighted";
				ar << unweighted;
			}
			ar << nEdges << sampleSize;
			if(_isWeighted)
			{
				ar << weights;
				std::string endWeights = "end_weights";
				ar << endWeights;
			}
			ar << data;
			if(nStoredBits > 0)
			{
				ar << storedBits;
			}
			std::string endFile = "end_distributions";
			ar << endFile;
		}
		template<class Archive> void load(Archive & ar, const unsigned int version)
		{
			std::string typeString;
			ar >> typeString;
			if(typeString != "empiricalDistribution")
			{
				throw std::runtime_error("File did not start with correct type specifier");
			}
			std::string weightString;
			ar >> weightString;
			if(weightString == "weighted")
			{
				_isWeighted = true;
			}
			else if(weightString == "unweighted")
			{
				_isWeighted = false;
			}
			else
			{
				throw std::runtime_error("Distributions file must start with either the string 'weighted' or 'unweighted'");
			}
			ar >> nEdges >> sampleSize;
			std::size_t bitsRequired = nEdges * sampleSize;
			std::size_t wholeIntsRequired = bitsRequired / (sizeof(int)*8);
			if(_isWeighted) 
			{
				ar >> weights;
				if(weights.size() != sampleSize)
				{
					throw std::runtime_error("Wrong number of weights loaded");
				}
				std::string endWeights;
				ar >> endWeights;
				if(endWeights != "end_weights")
				{
					throw std::runtime_error("Weights sections must end with the string 'end_weights'");
				}
			}
			ar >> data;
			if(data.size() != wholeIntsRequired)
			{
				throw std::runtime_error("Loaded data vector had wrong size");
			}
			int remainingBitsToRead = (int)(bitsRequired - (wholeIntsRequired * sizeof(int)*8));
			if(remainingBitsToRead > 0)
			{
				ar >> storedBits;
				nStoredBits = remainingBitsToRead;
				storedBits <<= (8*sizeof(int) - remainingBitsToRead);
			}
			else storedBits = nStoredBits = 0;
			std::string endDistributions;
			ar >> endDistributions;
			if(endDistributions != "end_distributions")
			{
				throw std::runtime_error("Distributions file must end with the string 'end_distributions'");
			}
		}
		BOOST_SERIALIZATION_SPLIT_MEMBER()
		void internalAdd(const EdgeState* state);
		empiricalDistribution();
		std::size_t nEdges;
		std::size_t sampleSize;
		bool _isWeighted;
		std::vector<double> weights;
		std::vector<int> data;
		int nStoredBits;
		int storedBits;
	};
}
#endif
