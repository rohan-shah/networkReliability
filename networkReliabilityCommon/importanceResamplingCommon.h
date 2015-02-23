#ifndef IMPORTANCE_RESAMPLING_COMMON_HEADER_GUARD
#define IMPORTANCE_RESAMPLING_COMMON_HEADER_GUARD
#include <vector>
#include "Context.h"
#include <boost/random/mersenne_twister.hpp>
#include "subObs/withImportanceResampling.h"
#include "NetworkReliabilityObsTree.h"
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/sum.hpp>
#include <boost/accumulators/statistics/stats.hpp>
namespace networkReliability
{
	struct importanceResamplingInput
	{
	public:
		importanceResamplingInput(const Context& context);
		bool shouldOutputTree;
		const Context& context;
		std::vector<double> thresholds;
		std::size_t n;
	private:
		importanceResamplingInput();
	};
	struct importanceResamplingOutput
	{
	public:
		importanceResamplingOutput(std::vector<::networkReliability::subObs::withImportanceResampling>& observations, boost::mt19937& randomSource, const Context& context, const std::vector<double>& thresholds);
		std::vector<::networkReliability::subObs::withImportanceResampling>& observations;
		boost::mt19937& randomSource;
		bool zeroEstimate;
		NetworkReliabilityObsTree tree;
		std::vector<boost::accumulators::accumulator_set<mpfr_class, boost::accumulators::stats<boost::accumulators::tag::sum> > > probabilities;
	private:
		importanceResamplingOutput();
	};
	void doImportanceResampling(const importanceResamplingInput& input, importanceResamplingOutput& output);
}
#endif
