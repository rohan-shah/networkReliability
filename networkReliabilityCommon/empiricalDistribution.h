#ifndef EMPIRICAL_DISTRIBUTION_HEADER_GUARD
#define EMPIRICAL_DISTRIBUTION_HEADER_GUARD
#include <boost/noncopyable.hpp>
#include "context.h"
#include "edgeState.h"
#include <string>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/vector.hpp>
#include <stdexcept>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "binaryDataSet.h"
namespace networkReliability
{
	class empiricalDistribution : protected binaryDataSet1, public boost::noncopyable
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
		empiricalDistribution(bool isWeighted, std::size_t nEdges, const context& context);
		void hintDataCount(std::size_t size);
		void add(const edgeState* state);
		void add(const edgeState* state, double weight);
		void expand(std::size_t count, std::vector<int>& output) const;
		std::size_t getNSamples() const;
		std::size_t getNEdges() const;
		bool isWeighted() const;
		double getWeight(std::size_t index) const;
		const context& getContext() const;
	private:
		template<class Archive> void save(Archive & ar, const unsigned int version) const
		{
			std::string typeString = "empiricalDistribution";
			ar << typeString;
			if(containedContext)
			{
				ar << *containedContext.get();
			}
			else ar << *externalContext;
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
			ar << *static_cast<const binaryDataSet1*>(this);
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
			containedContext.reset(new context(ar));
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
			ar >> *static_cast<binaryDataSet1*>(this);
			std::size_t bitsRequired = nEdges * sampleSize;
			std::size_t wholeIntsRequired = bitsRequired / (sizeof(int)*8);
			if(data.size() != wholeIntsRequired)
			{
				throw std::runtime_error("Loaded data vector had wrong size");
			}
			int remainingBitsToRead = (int)(bitsRequired - (wholeIntsRequired * sizeof(int)*8));
			int expectedNStoredBits;
			if(remainingBitsToRead > 0)
			{
				expectedNStoredBits = remainingBitsToRead;
			}
			else expectedNStoredBits = 0;
			if(nStoredBits != expectedNStoredBits)
			{
				throw std::runtime_error("Remaining number of stored bits had an unexpected value");
			}
			std::string endDistributions;
			ar >> endDistributions;
			if(endDistributions != "end_distributions")
			{
				throw std::runtime_error("Distributions file must end with the string 'end_distributions'");
			}
		}
		BOOST_SERIALIZATION_SPLIT_MEMBER()
		empiricalDistribution();
		std::size_t nEdges;
		std::size_t sampleSize;
		bool _isWeighted;
		std::vector<double> weights;

		std::shared_ptr<context> containedContext;
		context const* externalContext;
	};
}
#endif
