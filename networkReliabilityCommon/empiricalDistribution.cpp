#include "empiricalDistribution.h"
#include <stdexcept>
#include <fstream>
namespace networkReliability
{
	empiricalDistribution::empiricalDistribution(bool _isWeighted, std::size_t nEdges)
	:_isWeighted(_isWeighted), nStoredBits(0), storedBits(0), nEdges(nEdges), sampleSize(0)
	{}
	void empiricalDistribution::hintDataCount(std::size_t size)
	{
		if(_isWeighted) weights.reserve(size);
		std::size_t nBitsRequired = size*nEdges;
		std::size_t nWholeIntsRequired = nBitsRequired / (8*sizeof(int));
		if(nWholeIntsRequired * sizeof(int)*8 == nBitsRequired)
		{
			data.reserve(nWholeIntsRequired);
		}
		else data.reserve(nWholeIntsRequired+1);
	}
	void empiricalDistribution::internalAdd(const EdgeState* state)
	{
		for(int k = 0; k < nEdges; k++)
		{
			nStoredBits++;
			storedBits<<=1;
			storedBits += ((state[k] & OP_MASK) > 0);
			if(nStoredBits == sizeof(int)*8)
			{
				data.push_back(storedBits);
				storedBits = nStoredBits = 0;
			}
		}
		sampleSize++;
	}
	bool empiricalDistribution::isWeighted() const
	{
		return _isWeighted;
	}
	double empiricalDistribution::getWeight(std::size_t index) const
	{
		if(!_isWeighted) throw std::runtime_error("Trying to get weight of unweighted empirical distribution");
		return weights[index];
	}
	void empiricalDistribution::add(const EdgeState* state)
	{
		if(_isWeighted)
		{
			throw std::runtime_error("Attempting to add unweighted point to weighted empirical distribution");
		}
		internalAdd(state);
	}
	void empiricalDistribution::add(const EdgeState* state, double weight)
	{
		if(!_isWeighted)
		{
			throw std::runtime_error("Attempting to add weighted point to unweighted empirical distribution");
		}
		internalAdd(state);
		weights.push_back(weight);
	}
	empiricalDistribution& empiricalDistribution::operator=(empiricalDistribution&& other)
	{
		_isWeighted = other._isWeighted;
		weights.swap(other.weights);
		data.swap(other.data);
		nStoredBits = other.nStoredBits;
		storedBits = other.storedBits;
		sampleSize = other.sampleSize;
		return *this;
	}
	empiricalDistribution::empiricalDistribution(empiricalDistribution&& other)
	{
		_isWeighted = other._isWeighted;
		weights.swap(other.weights);
		data.swap(other.data);
		nStoredBits = other.nStoredBits;
		storedBits = other.storedBits;
		sampleSize = other.sampleSize;
	}
	void empiricalDistribution::expand(int count, std::vector<int>& output)
	{
		if(output.size() != nEdges) throw std::runtime_error("Wrong number of elements in vector passed to empiricalDistribution::expand");
		int initialBit = (int)(count*nEdges);
		int initialInt = (int)(initialBit / (sizeof(int)*8));
		int storedBits = data[initialInt];
		int extraBitsRead = initialBit - (initialInt* sizeof(int)*8);
		int nStoredBits = sizeof(int)*8;
		if(extraBitsRead > 0)
		{
			storedBits <<= extraBitsRead; 
			nStoredBits = sizeof(int)*8 - extraBitsRead;
		}
		int currentInt = initialInt;
		for(std::size_t edgeCounter = 0; edgeCounter < nEdges; edgeCounter++)
		{
			if(nStoredBits == 0)
			{
				currentInt++;
				storedBits = data[currentInt];
				nStoredBits = sizeof(int)*8;
			}
			if((storedBits & (1 << (sizeof(int)*8-1))) != 0)
			{
				output[edgeCounter] = 1;
			}
			else output[edgeCounter] = 0;
			storedBits<<=1;
			nStoredBits--;
		}
		
	}
	empiricalDistribution::empiricalDistribution()
	{}
	std::size_t empiricalDistribution::getNSamples() const
	{
		return sampleSize;
	}
	std::size_t empiricalDistribution::getNEdges() const
	{
		return nEdges;
	}

}
