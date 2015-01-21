#include "binaryDataSet.h"
#include <stdexcept>
namespace networkReliability
{
	void binaryDataSet1::add(const EdgeState* state, const std::size_t size)
	{
		for(int k = 0; k < size; k++)
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
		}
		throw std::runtime_error("Internal error");
	}
	void binaryDataSet2::add(const EdgeState* state, const std::size_t size)
	{
		for(int k = 0; k < size; k++)
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
	void binaryDataSet1::expand(int index, std::vector<int>& output) const 
	{
		const std::size_t nEdges = output.size();
		int initialBit = (int)(index*nEdges);
		int initialInt = (int)(initialBit / (sizeof(unsigned int)*8));
		unsigned int storedBits = data[initialInt];
		int extraBitsRead = initialBit - (initialInt* sizeof(unsigned int)*8);
		int nStoredBits = sizeof(unsigned int)*8;
		if(extraBitsRead > 0)
		{
			storedBits <<= extraBitsRead;
			nStoredBits = sizeof(unsigned int)*8 - extraBitsRead;
		}
		int currentInt = initialInt;
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
	void binaryDataSet2::expand(int index, EdgeState* output, const std::size_t nEdges) const
	{
		int initialBit = (int)(index*nEdges*2);
		int initialInt = (int)(initialBit / (sizeof(unsigned int)*8));
		unsigned int storedBits = data[initialInt];
		int extraBitsRead = initialBit - (initialInt* sizeof(unsigned int)*8);
		int nStoredBits = sizeof(unsigned int)*8;
		if(extraBitsRead > 0)
		{
			storedBits <<= extraBitsRead;
			nStoredBits = sizeof(unsigned int)*8 - extraBitsRead;
		}
		int currentInt = initialInt;
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
}
