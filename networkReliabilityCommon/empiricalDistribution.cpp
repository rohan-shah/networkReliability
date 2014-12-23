#include "empiricalDistribution.h"
#include <stdexcept>
#include <fstream>
namespace networkReliability
{
	empiricalDistribution::empiricalDistribution(bool isWeighted, std::size_t nEdges)
	:isWeighted(isWeighted), nStoredBits(0), storedBits(0), nEdges(nEdges), sampleSize(0)
	{}
	void empiricalDistribution::hintDataCount(std::size_t size)
	{
		if(isWeighted) weights.reserve(size);
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
	void empiricalDistribution::add(const EdgeState* state)
	{
		if(isWeighted)
		{
			throw std::runtime_error("Attempting to add unweighted point to weighted empirical distribution");
		}
		internalAdd(state);
	}
	void empiricalDistribution::add(const EdgeState* state, double weight)
	{
		if(!isWeighted)
		{
			throw std::runtime_error("Attempting to add weighted point to unweighted empirical distribution");
		}
		internalAdd(state);
		weights.push_back(weight);
	}
	empiricalDistribution& empiricalDistribution::operator=(empiricalDistribution&& other)
	{
		isWeighted = other.isWeighted;
		weights.swap(other.weights);
		data.swap(other.data);
		nStoredBits = other.nStoredBits;
		storedBits = other.storedBits;
		sampleSize = other.sampleSize;
	}
	empiricalDistribution::empiricalDistribution(empiricalDistribution&& other)
	{
		isWeighted = other.isWeighted;
		weights.swap(other.weights);
		data.swap(other.data);
		nStoredBits = other.nStoredBits;
		storedBits = other.storedBits;
		sampleSize = other.sampleSize;
	}
	empiricalDistribution empiricalDistribution::load(std::string file)
	{
		std::ifstream stream(file.c_str(), std::ios_base::binary);
		if(!stream)
		{
			throw std::runtime_error("Unable to load specified file");
		}
		empiricalDistribution loaded;
		std::string weightString;
		std::getline(stream, weightString, '\0');
		if(weightString == "weighted")
		{
			loaded.isWeighted = true;
		}
		else if(weightString == "unweighted")
		{
			loaded.isWeighted = false;
		}
		else
		{
			throw std::runtime_error("Distributions file must start with either the string 'weighted' or 'unweighted'");
		}
		stream.read((char*)&loaded.nEdges, sizeof(std::size_t));
		stream.read((char*)&loaded.sampleSize, sizeof(std::size_t));
		std::size_t bitsRequired = loaded.nEdges * loaded.sampleSize;
		std::size_t wholeIntsRequired = bitsRequired / (sizeof(int)*8);
		if(wholeIntsRequired * sizeof(int) * 8 == bitsRequired)
		{
			loaded.data.resize(wholeIntsRequired);
		}
		else 
		{
			loaded.data.resize(wholeIntsRequired+1);
		}
		if(loaded.isWeighted) 
		{
			loaded.weights.resize(loaded.sampleSize);
			stream.read((char*)&(loaded.weights[0]), sizeof(double)*loaded.sampleSize);
			std::string endWeights;
			stream >> endWeights;
			if(endWeights != "end_weights")
			{
				throw std::runtime_error("Weights sections must end with the string 'end_weights'");
			}
		}
		stream.read((char*)&(loaded.data[0]), sizeof(int)*wholeIntsRequired);
		int remainingBitsToRead = bitsRequired - (wholeIntsRequired * sizeof(int)*8);
		if(remainingBitsToRead > 0)
		{
			stream.read((char*)&loaded.storedBits, sizeof(int));
			loaded.nStoredBits = remainingBitsToRead;
			loaded.storedBits <<= (8*sizeof(int) - remainingBitsToRead);
		}
		std::string endDistributions;
		stream >> endDistributions;
		if(endDistributions != "end_distributions")
		{
			throw std::runtime_error("Distributions file must end with the string 'end_distributions'");
		}
		int position = stream.tellg();
		stream.seekg(0, stream.end);
		int length = stream.tellg();
		if(position != length)
		{
			throw std::runtime_error("Distributions file was too long");
		}
		return loaded;
	}
	void empiricalDistribution::expand(int count, std::vector<int>& output)
	{
		if(output.size() != nEdges) throw std::runtime_error("Wrong number of elements in vector passed to empiricalDistribution::expand");
		int initialBit = count*nEdges;
		int initialInt = initialBit / (sizeof(int)*8);
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
			if((storedBits & (1 << (sizeof(int)*8)-1)) != 0)
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
	bool empiricalDistribution::save(std::string file)
	{
		std::ofstream stream(file.c_str(), std::ios_base::binary);
		if(!stream)
		{
			return false;
		}
		if(isWeighted)
		{
			std::string weighted = "weighted";
			stream << weighted << '\0';
		}
		else
		{
			std::string unweighted = "unweighted";
			stream << unweighted << '\0';
		}
		stream.write((char*)&nEdges, sizeof(std::size_t));
		stream.write((char*)&sampleSize, sizeof(std::size_t));
		if(isWeighted)
		{
			stream.write((char*)&(weights[0]), sizeof(double)*weights.size());
			std::string endWeights = "end_weights";
			stream << endWeights;
		}
		stream.write((char*)&(data[0]), sizeof(int)*data.size());
		if(nStoredBits > 0)
		{
			stream.write((char*)&storedBits, sizeof(int));
		}
		std::string endFile = "end_distributions";
		stream << endFile;
		stream.flush();
		stream.close();
		return true;
	}
	std::size_t empiricalDistribution::getNSamples() const
	{
		return sampleSize;
	}
	std::size_t empiricalDistribution::getNEdges() const
	{
		return nEdges;
	}

}
