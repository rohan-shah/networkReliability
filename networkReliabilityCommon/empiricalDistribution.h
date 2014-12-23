#ifndef EMPIRICAL_DISTRIBUTION_HEADER_GUARD
#define EMPIRICAL_DISTRIBUTION_HEADER_GUARD
#include <boost/noncopyable.hpp>
#include "astCode.hpp"
#include "EdgeState.h"
#include <string>
namespace networkReliability
{
	class empiricalDistribution : public boost::noncopyable
	{
	public:
		empiricalDistribution(empiricalDistribution&& other);
		empiricalDistribution& operator=(empiricalDistribution&& other);
		empiricalDistribution(bool isWeighted, std::size_t nEdges);
		void hintDataCount(std::size_t size);
		bool save(std::string file);
		static empiricalDistribution load(std::string file);
		void add(const EdgeState* state);
		void add(const EdgeState* state, double weight);
		void expand(int count, std::vector<int>& output);
		std::size_t getNSamples() const;
		std::size_t getNEdges() const;
	private:
		void internalAdd(const EdgeState* state);
		empiricalDistribution();
		std::size_t nEdges;
		std::size_t sampleSize;
		bool isWeighted;
		std::vector<double> weights;
		std::vector<int> data;
		int nStoredBits;
		int storedBits;
	};
}
#endif
