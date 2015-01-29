#ifndef RESAMPLING_COMMON_HEADER_GUARD
#define RESAMPLING_COMMON_HEADER_GUARD
#include <vector>
#include "Context.h"
#include <boost/random/mersenne_twister.hpp>
#include "NetworkReliabilitySubObs.h"
#include "NetworkReliabilitySubObsTree.h"
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/sum.hpp>
namespace networkReliability
{
	struct resamplingInput
	{
	public:
		resamplingInput(const Context& context, bool compensateInResampling);
		bool shouldOutputTree;
		const Context& context;
		std::vector<double> thresholds;
		std::size_t n;
		int finalSplittingStep;
		bool compensateInResampling;
	private:
		resamplingInput();
	};
	struct resamplingOutput
	{
	public:
		resamplingOutput(std::vector<NetworkReliabilitySubObs>& observations, boost::mt19937& randomSource, const Context& context, const std::vector<double>& thresholds);
		std::vector<NetworkReliabilitySubObs>& observations;
		boost::mt19937& randomSource;
		bool zeroEstimate;
		NetworkReliabilitySubObsTree tree;
		std::vector<boost::accumulators::accumulator_set<mpfr_class, boost::accumulators::stats<boost::accumulators::tag::sum> > > probabilities;
	private:
		resamplingOutput();
	};
	void doResampling(const resamplingInput& input, resamplingOutput& output);
}
#endif
