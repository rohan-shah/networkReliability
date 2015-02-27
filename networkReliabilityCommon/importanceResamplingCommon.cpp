#include "importanceResamplingCommon.h"
#include "NetworkReliabilityObs.h"
#include "aliasMethod.h"
#include <boost/math/special_functions/binomial.hpp>
#include "obs/withImportanceResampling.h"
namespace networkReliability
{
	importanceResamplingInput::importanceResamplingInput(networkReliability::Context const& context)
		:context(context)
	{}
	importanceResamplingOutput::importanceResamplingOutput(std::vector<::networkReliability::subObs::withImportanceResampling>& observations, boost::mt19937& randomSource, const Context& context, const std::vector<double>& thresholds)
		:observations(observations), randomSource(randomSource), tree(&context, thresholds)
	{
		boost::accumulators::accumulator_set<mpfr_class, boost::accumulators::stats<boost::accumulators::tag::sum> > zeroInitialisedAccumulator(boost::parameter::keyword<boost::accumulators::tag::sample>::get() = 0);
		probabilities.resize(thresholds.size(), zeroInitialisedAccumulator);
	}
	void doImportanceResampling(const importanceResamplingInput& input, importanceResamplingOutput& output)
	{
		std::vector<int> componentsVector;
		output.zeroEstimate = false;
		//The ith entry in this vector gives the index within the n subObservations of the current generation, of the ith potentially disconnected subObservation. We need this for the tree.
		std::vector<int> potentiallyDisconnectedIndices, nextPotentiallyDisconnectedIndices;
		potentiallyDisconnectedIndices.reserve(input.n);
		nextPotentiallyDisconnectedIndices.reserve(input.n);

		std::vector<::networkReliability::subObs::withImportanceResampling> nextStepObservations;

		//Used in the alias method application
		std::vector<double> resamplingProbabilities;
		std::vector<std::ptrdiff_t> aliasMethodTemporary1, aliasMethodTemporary2;
		std::vector<std::pair<double, std::ptrdiff_t> > aliasMethodTemporary3;
		//Construct another version of the thresholds input, with an extra zero at the end
		std::vector<double> thresholds(input.thresholds.begin(), input.thresholds.end());
		thresholds.push_back(0);
		for (std::size_t i = 0; i < input.n; i++)
		{
			::networkReliability::obs::withImportanceResampling currentObs = ::networkReliability::obs::withImportanceResampling::constructConditional(input.context, output.randomSource);
			::networkReliability::subObs::withImportanceResampling subObs = ::networkReliability::obs::getSubObservation<::networkReliability::obs::withImportanceResampling>::get(currentObs, thresholds[0], thresholds[1]);

			if(input.shouldOutputTree) output.tree.add(subObs, 0, -1, subObs.getMinCut() < HIGH_CAPACITY);
			if (subObs.getMinCut() >= HIGH_CAPACITY)
			{
				output.probabilities[0](0);
			}
			else
			{
				output.probabilities[0](currentObs.getConditioningProb());
				nextStepObservations.push_back(std::move(subObs));
				nextPotentiallyDisconnectedIndices.push_back(i);
			}
		}
		if(nextStepObservations.size() == 0)
		{
			output.zeroEstimate = true;
			return;
		}
		for(std::size_t splittingLevel = 0; splittingLevel < input.thresholds.size()-1; splittingLevel++)
		{
			//resampling step
			output.observations.clear();
			potentiallyDisconnectedIndices.clear();
			resamplingProbabilities.clear();
			mpfr_class sum = 0;
			for (std::vector<::networkReliability::subObs::withImportanceResampling>::iterator j = nextStepObservations.begin(); j != nextStepObservations.end(); j++)
			{
				sum += j->getGeneratedObservationConditioningProb();
				resamplingProbabilities.push_back(j->getGeneratedObservationConditioningProb().convert_to<double>());
			}
			if(sum.convert_to<double>() == 0) throw std::runtime_error("Sum of importance weights was zero");
			mpfr_class averageWeight = sum / input.n;
			aliasMethod::aliasMethod alias(resamplingProbabilities, sum.convert_to<double>(), aliasMethodTemporary1, aliasMethodTemporary2, aliasMethodTemporary3);
			for (std::size_t k = 0; k < input.n; k++)
			{
				int index = alias(output.randomSource);
				output.observations.push_back(nextStepObservations[index].copyWithGeneratedObservationConditioningProb(averageWeight));
				potentiallyDisconnectedIndices.push_back(nextPotentiallyDisconnectedIndices[index]);
			}
			nextPotentiallyDisconnectedIndices.clear();
			nextStepObservations.clear();
			for(std::vector<::networkReliability::subObs::withImportanceResampling>::iterator j = output.observations.begin(); j != output.observations.end(); j++)
			{
				::networkReliability::obs::withImportanceResampling newObs = ::networkReliability::subObs::getObservation<::networkReliability::subObs::withImportanceResampling>::get(*j, output.randomSource);
				::networkReliability::subObs::withImportanceResampling sub = ::networkReliability::obs::getSubObservation<::networkReliability::obs::withImportanceResampling>::get(newObs, thresholds[splittingLevel + 1], thresholds[splittingLevel + 2]);
				if(input.shouldOutputTree) output.tree.add(sub, splittingLevel+1, potentiallyDisconnectedIndices[std::distance(output.observations.begin(), j)], sub.getMinCut() < HIGH_CAPACITY);

				if(sub.getMinCut() < HIGH_CAPACITY)
				{
					nextStepObservations.push_back(std::move(sub));
					output.probabilities[splittingLevel+1](newObs.getConditioningProb());
					nextPotentiallyDisconnectedIndices.push_back(std::distance(output.observations.begin(), j));
				}
			}
			if(nextStepObservations.size() == 0)
			{
				output.zeroEstimate = true;
				return;
			}
		}
		output.observations.swap(nextStepObservations);
	}
}
