#ifndef RESAMPLING_COMMON_HEADER_GUARD
#define RESAMPLING_COMMON_HEADER_GUARD
#include <vector>
#include "context.h"
#include <boost/random/mersenne_twister.hpp>
#include "networkReliabilityObsTree.h"
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/sum.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include "subObs/withResampling.h"
namespace networkReliability
{
	struct resamplingInput
	{
	public:
		resamplingInput(const context& contextObj);
		bool shouldOutputTree;
		const context& contextObj;
		std::vector<double> thresholds;
		std::size_t n;
		int finalSplittingStep;
	private:
		resamplingInput();
	};
	struct resamplingOutput
	{
	public:
		resamplingOutput(std::vector<::networkReliability::subObs::withResampling>& observations, boost::mt19937& randomSource, const context& contextObj, const std::vector<double>& thresholds);
		std::vector<::networkReliability::subObs::withResampling>& observations;
		boost::mt19937& randomSource;
		bool zeroEstimate;
		NetworkReliabilityObsTree tree;
		std::vector<boost::accumulators::accumulator_set<mpfr_class, boost::accumulators::stats<boost::accumulators::tag::sum> > > probabilities;
	private:
		resamplingOutput();
	};
	void doResampling(const resamplingInput& input, resamplingOutput& output);
}
#endif
