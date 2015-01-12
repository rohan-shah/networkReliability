#include "binaryDataSet.h"
namespace networkReliability
{
	void binaryDataSet1::add(const EdgeState* state, const std::size_t size)
	{
		for(int k = 0; k < size; k++)
		{
			nStoredBits++;
			storedBits <<=1;
			storedBits += ((state[k] & OP_MASK) > 0);
			if(nStoredBits == sizeof(int)*8)
			{
				data.push_back(storedBits);
				storedBits = nStoredBits = 0;
			}
		}
	}
	void binaryDataSet1::expand(int index, std::vector<int>& output)
	{
		const std::size_t nEdges = output.size();
		int initialBit = (int)(index*nEdges);
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
				if(currentInt == data.size())
				{
					storedBits = this->storedBits;
					nStoredBits = this->nStoredBits;
					storedBits <<= (8*sizeof(int) - nStoredBits);
				}
				else
				{
					storedBits = data[currentInt];
					nStoredBits = sizeof(int)*8;
				}
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
}
