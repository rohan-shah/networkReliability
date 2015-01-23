#include <boost/program_options.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include "Arguments.h"
#include "NetworkReliabilityObs.h"
#include "NetworkReliabilitySubObs.h"
#include "NetworkReliabilitySubObsCollection.h"
#include <vector>
#include "graphAlgorithms.h"
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/sum.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include "includeMPFR.h"
#include "ArgumentsMPFR.h"
#include "conditionalPMC.h"
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/special_functions.hpp>
#include <fstream>
#include <boost/random/uniform_real_distribution.hpp>
#include "aliasMethod.h"
#include "empiricalDistribution.h"
#include "depth_first_search_restricted.hpp"
#include "connected_components_restricted.hpp"
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/math/distributions.hpp>
#include "NetworkReliabilitySubObsTree.h"
namespace networkReliability
{
	int main(int argc, char **argv)
	{
		typedef NetworkReliabilityObs::conditioning_type calculation_type;

		boost::program_options::options_description options("Usage");
		options.add_options()
			("gridGraph", boost::program_options::value<int>(), "(int) The dimension of the square grid graph to use. Incompatible with graphFile. ")
			("graphFile", boost::program_options::value<std::string>(), "(string) The path to a graphml file. Incompatible with gridGraph")
			("opProbability", boost::program_options::value<std::string>(), "(float) The probability that an edge is operational")
			("seed", boost::program_options::value<int>(), "(int) The random seed used to generate the random graphs")
			("interestVertices", boost::program_options::value<std::vector<int> >()->multitoken(), "(int) The vertices of interest, that should be connected. ")
			("initialRadius", boost::program_options::value<int>(), "(int) The initial radius to use")
			("n", boost::program_options::value<std::size_t>(), "(int) The number of graphs initially generated")
			("outputConditionalDistribution", boost::program_options::value<std::string>(), "(path) File to output the empirical conditional distribution")
			("outputTree", boost::program_options::value<std::string>(), "(path) File to output simulation tree to")
			("useSpatialDistances", boost::program_options::value<std::vector<double> >()->multitoken(), "(float) Input spatial distances must consist of two numbers; A maximum distance and the number of steps to take.")
			("help", "Display this message");

		boost::program_options::variables_map variableMap;
		try
		{
			boost::program_options::store(boost::program_options::parse_command_line(argc, argv, options), variableMap);
		}
		catch(boost::program_options::error& ee)
		{
			std::cerr << "Error parsing command line arguments: " << ee.what() << std::endl << std::endl;
			std::cerr << options << std::endl;
			return -1;
		}
		if(variableMap.count("help") > 0)
		{
			std::cout << options << std::endl;
			return 0;
		}

		mpfr_class opProbability;
		if(!readProbabilityString(variableMap, opProbability))
		{
			std::cout << "Unable to read input opProbability" << std::endl;
			return 0;
		}
		mpfr_class inopProbability = 1 - opProbability;

		Context context = Context::emptyContext();
		if(!readContext(variableMap, context, opProbability))
		{
			std::cout << "Unable to construct context object" << std::endl;
			return 0;
		}
		context.setMinCut(true);
		const std::size_t nEdges = context.getNEdges();

		std::vector<double> thresholds;
		std::string message;
		if(!readThresholds(variableMap, thresholds, message))
		{
			std::cout << message << std::endl;
			return 0;
		}

		std::size_t n;
		if(!readN(variableMap, n))
		{
			return 0;
		}

		boost::mt19937 randomSource;
		if(variableMap.count("seed") > 0)
		{
			randomSource.seed(variableMap["seed"].as<int>());
		}

		std::vector<NetworkReliabilitySubObs> observations;
		std::vector<NetworkReliabilitySubObs> nextStepObservations;

		boost::accumulators::accumulator_set<calculation_type, boost::accumulators::stats<boost::accumulators::tag::sum> > zeroInitialisedAccumulator(boost::parameter::keyword<boost::accumulators::tag::sample>::get() = 0);
		//mean of getting to the next level, once we've conditioned on having enough edges.
		std::vector<boost::accumulators::accumulator_set<calculation_type, boost::accumulators::stats<boost::accumulators::tag::sum> > > probabilities(thresholds.size(), zeroInitialisedAccumulator);

		//working data for graph algorithms
		std::vector<int> components;
		std::vector<boost::default_color_type> colorMap;
		boost::detail::depth_first_visit_restricted_impl_helper<Context::internalGraph>::stackType stack;

		//Working data for alias method
		std::vector<std::ptrdiff_t> aliasMethodTemporary1, aliasMethodTemporary2;
		std::vector<std::pair<double, std::ptrdiff_t> > aliasMethodTemporary3;

		calculation_type estimate = 0;

		int finalSplittingStep = thresholds.size()-1;

		//If we specify the useMinCut option then we need to do resampling. This vector holds the probabilities
		std::vector<double> resamplingProbabilities;

		//Object to hold the tree of generated objects, if input outputTree is specified
		bool outputTree = variableMap.count("outputTree") > 0;
		NetworkReliabilitySubObsTree tree(&context, thresholds.size(), thresholds);
		for (std::size_t i = 0; i < n; i++)
		{
			NetworkReliabilityObs currentObs = NetworkReliabilityObs::constructConditional(context, randomSource);
			NetworkReliabilitySubObs subObs = currentObs.getSubObservation(thresholds[0]);

			if (subObs.getMinCut() >= HIGH_CAPACITY)
			{
				probabilities[0](0);
				if(outputTree) tree.add(subObs, 0, -1, false);
			}
			else
			{
				probabilities[0](1);
				observations.push_back(std::move(subObs));
				if(outputTree) tree.add(subObs, 0, -1, true);
			}
		}
		if(observations.size() == 0)
		{
			estimate = 0;
			goto returnEstimate;
		}
		for(int splittingLevel = 0; splittingLevel < finalSplittingStep; splittingLevel++)
		{
			nextStepObservations.clear();
			for(std::vector<NetworkReliabilitySubObs>::iterator j = observations.begin(); j != observations.end(); j++)
			{
				NetworkReliabilityObs newObs = j->getObservation(randomSource);
				NetworkReliabilitySubObs sub = newObs.getSubObservation(thresholds[splittingLevel + 1]);
				if(sub.getMinCut() < HIGH_CAPACITY)
				{
					probabilities[splittingLevel+1](newObs.getConditioningProb());
					nextStepObservations.push_back(std::move(sub));
				}
			}
			if(nextStepObservations.size() == 0)
			{
				estimate = 0;
				goto returnEstimate;
			}
			std::size_t previousSize = observations.size();
			//resampling step. Don't do this on the last step, because we're already at radius 0
			if(splittingLevel < finalSplittingStep-1)
			{
				observations.clear();
				resamplingProbabilities.clear();
				mpfr_class sum = 0;
				for (std::vector<NetworkReliabilitySubObs>::iterator j = nextStepObservations.begin(); j != nextStepObservations.end(); j++)
				{
					sum += j->getGeneratedObservationConditioningProb();
					resamplingProbabilities.push_back(j->getGeneratedObservationConditioningProb().convert_to<double>());
				}
				if(sum.convert_to<double>() == 0) throw std::runtime_error("Sum of importance weights was zero");
				mpfr_class averageWeight = sum / previousSize;
				aliasMethod::aliasMethod alias(resamplingProbabilities, sum.convert_to<double>(), aliasMethodTemporary1, aliasMethodTemporary2, aliasMethodTemporary3);
				for (std::size_t k = 0; k < previousSize; k++)
				{
					int index = alias(randomSource);
					observations.push_back(nextStepObservations[index].copyWithGeneratedObservationConditioningProb(averageWeight));
				}
			}
			else observations.swap(nextStepObservations);
		}
		estimate = (boost::accumulators::sum(probabilities[finalSplittingStep]) / n);
		if (variableMap.count("outputConditionalDistribution") > 0)
		{
			empiricalDistribution outputDistributions(true, nEdges);
			for(std::vector<NetworkReliabilitySubObs>::iterator i = observations.begin(); i != observations.end(); i++)
			{
				outputDistributions.add(i->getState(), i->getConditioningProb().convert_to<double>());
			}
			try
			{
				std::ofstream outputStream(variableMap["outputConditionalDistribution"].as<std::string>().c_str(), std::ios_base::binary);
				boost::archive::binary_oarchive outputArchive(outputStream, boost::archive::no_codecvt);
				outputArchive << outputDistributions;
			}
			catch(std::runtime_error& err)
			{
				std::cout << "Error saving empirical distributions to file: " << err.what() << std::endl;
				return 0;
			}
		}
returnEstimate:
		std::cout << "Estimate is " << estimate.convert_to<double>() << std::endl;
		return 0;
	}
}
int main(int argc, char **argv)
{
	return networkReliability::main(argc, argv);
}
