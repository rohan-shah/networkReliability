#include <boost/program_options.hpp>
#include "Arguments.h"
#include "NetworkReliabilityObs.h"
#include "NetworkReliabilityObsCollection.h"
#include "subObs/withImportanceResampling.h"
#include "obs/withImportanceResampling.h"
#include <vector>
#include "graphAlgorithms.h"
#include <boost/random/bernoulli_distribution.hpp>
#include "includeMPFR.h"
#include "ArgumentsMPFR.h"
#include "conditionalPMC.h"
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/special_functions.hpp>
#include <fstream>
#include <boost/random/uniform_real_distribution.hpp>
#include "importanceResamplingCommon.h"
#include "empiricalDistribution.h"
#include "depth_first_search_restricted.hpp"
#include "connected_components_restricted.hpp"
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/math/distributions.hpp>
#include "NetworkReliabilityObsTree.h"
#include "commonOptions.h"
namespace networkReliability
{
	int main(int argc, char **argv)
	{
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
			("useSpatialDistances", boost::program_options::value<std::vector<double> >()->multitoken(), "(float) Input spatial distances must consist of two or three numbers numbers; A maximum distance, an optional minimum distance and the number of steps to take.")
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
		const std::vector<int>& interestVertices = context.getInterestVertices();
		if(interestVertices.size() > 2) 
		{
			std::cout << "Importance resampling can only be used for two-terminal network reliability" << std::endl;
			return 0;
		}
		const std::size_t nEdges = context.getNEdges();

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

		importanceResamplingInput resamplingInputs(context);
		resamplingInputs.shouldOutputTree = variableMap.count("outputTree") > 0;
		resamplingInputs.n = n;

		std::string message;
		if(!readThresholds(variableMap, resamplingInputs.thresholds, message))
		{
			std::cout << message << std::endl;
			return 0;
		}
		std::vector<::networkReliability::subObs::withImportanceResampling> observations;

		//working data for graph algorithms
		std::vector<int> components;
		std::vector<boost::default_color_type> colorMap;
		boost::detail::depth_first_visit_restricted_impl_helper<Context::internalGraph>::stackType stack;

		mpfr_class estimate = 0;

		importanceResamplingOutput resamplingOutputs(observations, randomSource, context, resamplingInputs.thresholds);
		doImportanceResampling(resamplingInputs, resamplingOutputs);
		if(resamplingOutputs.zeroEstimate)
		{
			estimate = 0;
			goto returnEstimate;
		}
		estimate = (boost::accumulators::sum(resamplingOutputs.probabilities[resamplingInputs.thresholds.size()-1]) / n);
		if (variableMap.count("outputConditionalDistribution") > 0)
		{
			empiricalDistribution outputDistributions(true, nEdges, context);
			for(std::vector<::networkReliability::subObs::withImportanceResampling>::iterator i = observations.begin(); i != observations.end(); i++)
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
		if(resamplingInputs.shouldOutputTree)
		{
			std::cout << "Beginning tree layout....";
			std::cout.flush();
			resamplingOutputs.tree.layout();
			std::cout << "Done" << std::endl;
			try
			{
				std::ofstream outputStream(variableMap["outputTree"].as<std::string>().c_str(), std::ios_base::binary);
				boost::archive::binary_oarchive outputArchive(outputStream, boost::archive::no_codecvt);
				outputArchive << resamplingOutputs.tree;
			}
			catch(std::runtime_error& err)
			{
				std::cout << "Error saving tree to file: " << err.what() << std::endl;
				return 0;
			}
		}
		std::cout << "Estimate is " << estimate.convert_to<double>() << std::endl;
		return 0;
	}
}
int main(int argc, char **argv)
{
	return networkReliability::main(argc, argv);
}
