#include "resamplingCommon.h"
#include "aliasMethod.h"
#include <boost/math/special_functions/binomial.hpp>
#include "obs/withResampling.h"
#include "subObs/withResampling.h"
#include "subObs/getObservation.hpp"
#include "obs/getSubObservation.hpp"
namespace networkReliability
{
	resamplingInput::resamplingInput(networkReliability::Context const& context, bool compensateInResampling)
		:context(context), compensateInResampling(compensateInResampling)
	{}
	resamplingOutput::resamplingOutput(std::vector<::networkReliability::subObs::withResampling>& observations, boost::mt19937& randomSource, const Context& context, const std::vector<double>& thresholds)
		:observations(observations), randomSource(randomSource), tree(&context, thresholds)
	{
		boost::accumulators::accumulator_set<mpfr_class, boost::accumulators::stats<boost::accumulators::tag::sum> > zeroInitialisedAccumulator(boost::parameter::keyword<boost::accumulators::tag::sample>::get() = 0);
		probabilities.resize(thresholds.size(), zeroInitialisedAccumulator);
	}
	void doResampling(const resamplingInput& input, resamplingOutput& output)
	{
		output.zeroEstimate = false;
		//The ith entry in this vector gives the index within the n subObservations of the current generation, of the ith potentially disconnected subObservation. We need this for the tree.
		std::vector<int> potentiallyDisconnectedIndices, nextPotentiallyDisconnectedIndices;
		potentiallyDisconnectedIndices.reserve(input.n);
		nextPotentiallyDisconnectedIndices.reserve(input.n);

		std::vector<::networkReliability::subObs::withResampling> nextStepObservations;

		//Used in the alias method application
		std::vector<double> resamplingProbabilities;
		std::vector<std::ptrdiff_t> aliasMethodTemporary1, aliasMethodTemporary2;
		std::vector<std::pair<double, std::ptrdiff_t> > aliasMethodTemporary3;

		for (std::size_t i = 0; i < input.n; i++)
		{
			::networkReliability::obs::withResampling currentObs = ::networkReliability::obs::withResampling::constructConditional(input.context, output.randomSource);
			::networkReliability::subObs::withResampling subObs = ::networkReliability::obs::getSubObservation<::networkReliability::obs::withResampling>::get(currentObs, input.thresholds[0]);

			if(input.shouldOutputTree) output.tree.add(subObs, 0, -1, subObs.getMinCut() < HIGH_CAPACITY);
			if (subObs.getMinCut() >= HIGH_CAPACITY)
			{
				output.probabilities[0](0);
			}
			else
			{
				output.probabilities[0](1);
				nextStepObservations.push_back(std::move(subObs));
				nextPotentiallyDisconnectedIndices.push_back((int)i);
			}
		}
		if(nextStepObservations.size() == 0)
		{
			output.zeroEstimate = true;
			return;
		}
		for(int splittingLevel = 0; splittingLevel < input.finalSplittingStep; splittingLevel++)
		{
			//resampling step
			output.observations.clear();
			potentiallyDisconnectedIndices.clear();
			resamplingProbabilities.clear();
			mpfr_class sum = 0;
			if(!input.compensateInResampling)
			{
				for (std::vector<::networkReliability::subObs::withResampling>::iterator j = nextStepObservations.begin(); j != nextStepObservations.end(); j++)
				{
					sum += j->getGeneratedObservationConditioningProb();
					resamplingProbabilities.push_back(j->getGeneratedObservationConditioningProb().convert_to<double>());
				}
			}
			else
			{
				for (std::vector<::networkReliability::subObs::withResampling>::iterator j = nextStepObservations.begin(); j != nextStepObservations.end(); j++)
				{
					int nUnfixed = (int)j->getPotentiallyDeactivated().size();
					int minCut = (int)j->getMinCut();
					sum += j->getGeneratedObservationConditioningProb() / boost::math::binomial_coefficient<mpfr_class>(nUnfixed, minCut);
					mpfr_class weight = j->getGeneratedObservationConditioningProb() / boost::math::binomial_coefficient<mpfr_class>(nUnfixed, minCut);
					resamplingProbabilities.push_back(weight.convert_to<double>());
				}
			}
			if(sum.convert_to<double>() == 0) throw std::runtime_error("Sum of importance weights was zero");
			mpfr_class averageWeight = sum / input.n;
			aliasMethod::aliasMethod alias(resamplingProbabilities, sum.convert_to<double>(), aliasMethodTemporary1, aliasMethodTemporary2, aliasMethodTemporary3);
			if(!input.compensateInResampling)
			{
				for (std::size_t k = 0; k < input.n; k++)
				{
					int index = (int)alias(output.randomSource);
					output.observations.push_back(nextStepObservations[index].copyWithGeneratedObservationConditioningProb(averageWeight));
					potentiallyDisconnectedIndices.push_back(nextPotentiallyDisconnectedIndices[index]);
				}
			}
			else
			{
				for (std::size_t k = 0; k < input.n; k++)
				{
					int index = (int)alias(output.randomSource);
					int nUnfixed = (int)nextStepObservations[index].getPotentiallyDeactivated().size();
					int minCut = (int)nextStepObservations[index].getMinCut();
					output.observations.push_back(nextStepObservations[index].copyWithGeneratedObservationConditioningProb(averageWeight * boost::math::binomial_coefficient<mpfr_class>(nUnfixed, minCut)));
					potentiallyDisconnectedIndices.push_back(nextPotentiallyDisconnectedIndices[index]);
				}
			}

			nextPotentiallyDisconnectedIndices.clear();
			nextStepObservations.clear();
			for(std::vector<::networkReliability::subObs::withResampling>::iterator j = output.observations.begin(); j != output.observations.end(); j++)
			{
				::networkReliability::obs::withResampling newObs = ::networkReliability::subObs::getObservation<::networkReliability::subObs::withResampling>::get(*j, output.randomSource);
				::networkReliability::subObs::withResampling sub = ::networkReliability::obs::getSubObservation<::networkReliability::obs::withResampling>::get(newObs, input.thresholds[splittingLevel + 1]);
				if(input.shouldOutputTree) output.tree.add(sub, splittingLevel+1, potentiallyDisconnectedIndices[std::distance(output.observations.begin(), j)], sub.getMinCut() < HIGH_CAPACITY);

				if(sub.getMinCut() < HIGH_CAPACITY)
				{
					nextStepObservations.push_back(std::move(sub));
					output.probabilities[splittingLevel+1](newObs.getConditioningProb());
					nextPotentiallyDisconnectedIndices.push_back((int)std::distance(output.observations.begin(), j));
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
