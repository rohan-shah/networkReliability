#include "binaryDataSet.h"
#include <stdexcept>
namespace networkReliability
{
	void binaryDataSet1::add(const EdgeState* state, const std::size_t size)
	{
		for(std::size_t k = 0; k < size; k++)
		{
			nStoredBits++;
			storedBits <<=1;
			storedBits += (unsigned int)((state[k] & OP_MASK) > 0);
			if(nStoredBits == sizeof(unsigned int)*8)
			{
				data.push_back(storedBits);
				storedBits = nStoredBits = 0u;
			}
		}
	}
	void binaryDataSet1::reserve(std::size_t nBits)
	{
		data.reserve((nBits / sizeof(unsigned int))+1);
	}
	void binaryDataSet2::reserve(std::size_t nStates)
	{
		data.reserve((2*nStates / sizeof(unsigned int))+1);
	}
	binaryDataSet::binaryDataSet(binaryDataSet&& other)
	{
		storedBits = other.storedBits;
		nStoredBits = other.nStoredBits;
		data.swap(other.data);
	}
	//A version more suitable for being saved to disk
	enum SavedEdgeState
	{
		FIXED_OP_SAVED = 0, FIXED_INOP_SAVED = 1, UNFIXED_OP_SAVED = 2, UNFIXED_INOP_SAVED = 3
	};
	inline EdgeState savedToStandard(SavedEdgeState saved)
	{
		switch(saved)
		{
		case FIXED_OP_SAVED:
			return FIXED_OP;
		case FIXED_INOP_SAVED:
			return FIXED_INOP;
		case UNFIXED_OP_SAVED:
			return UNFIXED_OP;
		case UNFIXED_INOP_SAVED:
			return UNFIXED_INOP;
		}
		throw std::runtime_error("Internal error");
	}
	inline SavedEdgeState standardToSaved(EdgeState state)
	{
		switch(state)
		{
		case FIXED_OP:
			return FIXED_OP_SAVED;
		case FIXED_INOP:
			return FIXED_INOP_SAVED;
		case UNFIXED_OP:
			return UNFIXED_OP_SAVED;
		case UNFIXED_INOP:
			return UNFIXED_INOP_SAVED;
		default:
			;
		}
		throw std::runtime_error("Internal error");
	}
	void binaryDataSet2::add(const EdgeState* state, const std::size_t size)
	{
		for(std::size_t k = 0; k < size; k++)
		{
			nStoredBits+=2;
			storedBits <<=2;
			storedBits += (unsigned int)standardToSaved(state[k]);
			if(nStoredBits == sizeof(unsigned int)*8)
			{
				data.push_back(storedBits);
				storedBits = nStoredBits = 0u;
			}
		}
	}
	void binaryDataSet1::expand(std::size_t index, std::vector<int>& output) const 
	{
		const std::size_t nEdges = output.size();
		std::size_t initialBit = index*nEdges;
		std::size_t initialInt = initialBit / (sizeof(unsigned int)*8);
		unsigned int storedBits = data[initialInt];
		std::size_t extraBitsRead = (std::size_t)(initialBit - (initialInt* sizeof(unsigned int)*8));
		unsigned int nStoredBits = (unsigned int)sizeof(unsigned int)*8;
		if(extraBitsRead > 0)
		{
			storedBits <<= extraBitsRead;
			nStoredBits = (unsigned int)(sizeof(unsigned int)*8 - extraBitsRead);
		}
		std::size_t currentInt = initialInt;
		for(std::size_t edgeCounter = 0; edgeCounter < nEdges; edgeCounter++)
		{
			if(nStoredBits == 0)
			{
				currentInt++;
				if(currentInt == data.size())
				{
					storedBits = this->storedBits;
					nStoredBits = this->nStoredBits;
					storedBits <<= (8*sizeof(unsigned int) - nStoredBits);
				}
				else
				{
					storedBits = data[currentInt];
					nStoredBits = sizeof(unsigned int)*8;
				}
			}
			if((storedBits & (1u << (sizeof(unsigned int)*8-1))) != 0)
			{
				output[edgeCounter] = 1;
			}
			else output[edgeCounter] = 0;
			storedBits<<=1;
			nStoredBits--;
		}
	}
	void binaryDataSet2::expand(std::size_t index, EdgeState* output, const std::size_t nEdges) const
	{
		std::size_t initialBit = index*nEdges*2;
		std::size_t initialInt = initialBit / (sizeof(unsigned int)*8);
		unsigned int storedBits = data[initialInt];
		std::size_t extraBitsRead = (std::size_t)(initialBit - (initialInt* sizeof(unsigned int)*8));
		std::size_t nStoredBits = sizeof(unsigned int)*8;
		if(extraBitsRead > 0)
		{
			storedBits <<= extraBitsRead;
			nStoredBits = sizeof(unsigned int)*8 - extraBitsRead;
		}
		std::size_t currentInt = initialInt;
		for(std::size_t edgeCounter = 0; edgeCounter < nEdges; edgeCounter++)
		{
			if(nStoredBits == 0)
			{
				currentInt++;
				if(currentInt == data.size())
				{
					storedBits = this->storedBits;
					nStoredBits = this->nStoredBits;
					storedBits <<= (8*sizeof(unsigned int) - nStoredBits);
				}
				else
				{
					storedBits = data[currentInt];
					nStoredBits = sizeof(unsigned int)*8;
				}
			}
			output[edgeCounter] = savedToStandard((SavedEdgeState)((storedBits & (3u << (sizeof(unsigned int)*8-2))) >> (sizeof(unsigned int)*8-2)));
			storedBits<<=2;
			nStoredBits-=2;
		}
	}
	binaryDataSet1& binaryDataSet1::operator=(binaryDataSet1&& other)
	{
		*static_cast<binaryDataSet*>(this) = static_cast<binaryDataSet&&>(other);
		return *this;
	}
	binaryDataSet& binaryDataSet::operator=(binaryDataSet&& other)
	{
		storedBits = other.storedBits;
		nStoredBits = other.nStoredBits;
		data.swap(other.data);
		return *this;
	}
	binaryDataSet1::binaryDataSet1(binaryDataSet1&& other)
		:binaryDataSet(std::move(other))
	{}
	binaryDataSet2::binaryDataSet2(binaryDataSet2&& other)
		:binaryDataSet(std::move(other))
	{}
	binaryDataSet2& binaryDataSet2::operator=(binaryDataSet2&& other)
	{
		*static_cast<binaryDataSet*>(this) = static_cast<binaryDataSet&&>(other);
		return *this;
	}
}
