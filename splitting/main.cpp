#include <boost/program_options.hpp>
#include "Arguments.h"
#include "NetworkReliabilityObs.h"
#include "obs/basic.h"
#include "subObs/basic.h"
#include <vector>
#include "graphAlgorithms.h"
#include <boost/accumulators/statistics/stats.hpp>
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
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
namespace networkReliability
{
	std::string toString(mpfr_class number)
	{
		mp_exp_t exponent;
		char* resultCStr = mpfr_get_str(NULL, &exponent, 10, 10, number.backend().data(), MPFR_RNDN);
		std::string resultStr  = resultCStr;
		free(resultCStr);
		return resultStr.substr(0, 1) + "." + resultStr.substr(1, resultStr.size() - 1) + "e" + boost::lexical_cast<std::string>(exponent - 1);
	}
	std::string toString(double number);
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
			("splittingFactor", boost::program_options::value<std::vector<float> >()->multitoken(), "(float) The splitting factor to use at every step")
			("n", boost::program_options::value<std::size_t>(), "(int) The number of graphs initially generated")
			("usePMC", "(Flag) Use PMC for the last step")
			("outputConditionalDistribution", boost::program_options::value<std::string>(), "(path) File to output the edge distribution of the conditional distribution")
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

		Context context = Context::emptyContext();
		if(!readContext(variableMap, context, opProbability))
		{
			std::cout << "Unable to construct context object" << std::endl;
			return 0;
		}

		bool usePMC = variableMap.count("usePMC") > 0;

		std::string message;
		std::vector<double> thresholds;
		if(!readThresholds(variableMap, thresholds, message))
		{
			std::cout << message << std::endl;
			return 0;
		}
		if(variableMap.count("splittingFactor") != 1)
		{
			std::cout << "Please enter a value for input splittingFactor" << std::endl;
			return 0;
		}

		std::vector<float> splittingFactors = variableMap["splittingFactor"].as<std::vector<float> >();
		if(splittingFactors.size() == 1)
		{
			splittingFactors.insert(splittingFactors.end(), thresholds.size() - 2, splittingFactors[0]);
		}
		if(splittingFactors.size() != thresholds.size() - 1)
		{
			std::cout << "Wrong number of values entered for input splittingFactor" << std::endl;
			return 0;
		}
		bool validSplittingFactors = true;

		for(std::vector<float>::const_iterator i = splittingFactors.begin(); i != splittingFactors.end(); i++)
		{
			if(*i < 1) validSplittingFactors = false;
		}
		if(!validSplittingFactors)
		{
			std::cout << "Input splittingFactors must contain numbers greater than 1" << std::endl;
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

		std::vector<::networkReliability::subObs::basic> observations;
		std::vector<::networkReliability::subObs::basic> nextStepObservations;

		const std::vector<int>& interestVertices = context.getInterestVertices();

		//working data for algorithms
		std::vector<int> components;
		std::vector<boost::default_color_type> colorMap;
		boost::detail::depth_first_visit_restricted_impl_helper<Context::internalGraph>::stackType stack;
		double totalSamples = (double)n;
		for (std::vector<float>::iterator k = splittingFactors.begin(); k != splittingFactors.end(); k++) totalSamples *= *k;
		std::size_t nEdges = context.getNEdges();

		mpfr_class firstMomentSum = 0;
		//If the usePMC flag is set, don't use splitting on the last step. Instead use PMC. 
		int finalSplittingThresholdIndex;
		if (usePMC) finalSplittingThresholdIndex = (int)(thresholds.size()-2);
		else finalSplittingThresholdIndex = (int)(thresholds.size()-1);

		//additional working data for getReducedGraph
		std::vector<int> edgeCounts;
		std::vector<int> reducedGraphInterestVertices(interestVertices.size());
		ConditionalTurnipInput turnipInput(randomSource, NULL, reducedGraphInterestVertices);
		turnipInput.exponentialRate = -boost::multiprecision::log(mpfr_class(1 - opProbability));

		//Work out the probability we initially condition on below by calling constructConditional
		const TruncatedBinomialDistribution::TruncatedBinomialDistribution& dist = context.getInopDistribution(0, nEdges, nEdges);
		const mpfr_class initialConditioningProb = 1 - dist.getCumulativeProbability((int)(context.getMinCutEdges() - 1));
		for (std::size_t i = 0; i < n; i++)
		{
			::networkReliability::obs::basic currentObs = ::networkReliability::obs::basic::constructConditional(context, randomSource);
			::networkReliability::subObs::basic subObs = ::networkReliability::obs::getSubObservation<::networkReliability::obs::basic>::get(currentObs, thresholds[0]);

			if (!isSingleComponent(context, subObs.getState(), components, stack, colorMap))
			{
				observations.push_back(std::move(subObs));
			}
		}
		for(int splittingLevel = 0; splittingLevel < finalSplittingThresholdIndex; splittingLevel++)
		{
			nextStepObservations.clear();
			boost::random::bernoulli_distribution<float> bernoulli(splittingFactors[splittingLevel] - floor(splittingFactors[splittingLevel]));
			for(std::vector<::networkReliability::subObs::basic>::iterator j = observations.begin(); j != observations.end(); j++)
			{
				int integerSplittingFactor = (int)floor(splittingFactors[splittingLevel]) + bernoulli(randomSource);
				for(int k = 0; k < integerSplittingFactor; k++)
				{
					::networkReliability::obs::basic newObs = ::networkReliability::subObs::getObservation<::networkReliability::subObs::basic>::get(*j, randomSource);
					::networkReliability::subObs::basic sub = ::networkReliability::obs::getSubObservation<::networkReliability::obs::basic>::get(newObs, thresholds[splittingLevel+1]);
					if(!isSingleComponent(context, sub.getState(), components, stack, colorMap))
					{
						nextStepObservations.push_back(std::move(sub));
					}
				}
			}
			observations.swap(nextStepObservations);
		}
		if (usePMC)
		{
			boost::random::bernoulli_distribution<float> bernoulli(*splittingFactors.rbegin() - floor(*splittingFactors.rbegin()));
			for (std::vector<::networkReliability::subObs::basic>::iterator j = observations.begin(); j != observations.end(); j++)
			{
				Context::internalGraph reducedGraph;
				j->getReducedGraph(reducedGraph, edgeCounts, components, stack, colorMap);
				turnipInput.n = (int)*splittingFactors.rbegin() + bernoulli(randomSource);
				turnipInput.graph = &reducedGraph;
				turnipInput.minimumInoperative = std::max(0, ((int)context.getMinCutEdges() - (int)j->getFixedInopCount()));
				const std::size_t nReducedVertices = boost::num_vertices(reducedGraph);
				const std::size_t nReducedEdges = boost::num_edges(reducedGraph);
				//If the graph is DEFINITELY disconnected, add a value of 1 and continue.
				if (nReducedEdges == 0)
				{
					bool alwaysConnected = true;
					for (std::size_t interestCounter = 1; interestCounter < interestVertices.size(); interestCounter++)
					{
						if (components[interestVertices[interestCounter]] != components[interestVertices[0]])
						{
							alwaysConnected = false;
							break;
						}
					}
					if (alwaysConnected) throw std::runtime_error("Internal error");
					firstMomentSum += turnipInput.n;
					continue;
				}
				turnipInput.edges.clear();
				turnipInput.edges.resize(nReducedEdges);
				turnipInput.edgeCounts.resize(nReducedEdges);
				Context::internalGraph::edge_iterator current, end;
				boost::tie(current, end) = boost::edges(reducedGraph);
				for (; current != end; current++)
				{
					int edgeIndex = boost::get(boost::edge_index, reducedGraph, *current);
					if (current->m_target != current->m_source)
					{
						turnipInput.edgeCounts[edgeIndex] = edgeCounts[current->m_source + current->m_target * nReducedVertices] + edgeCounts[current->m_target + current->m_source * nReducedVertices];
					}
					else
					{
						turnipInput.edgeCounts[edgeIndex] = edgeCounts[current->m_source + current->m_target * nReducedVertices];
					}
					turnipInput.edges[edgeIndex] = std::make_pair((int)current->m_source, (int)current->m_target);
				}
				for (std::size_t k = 0; k < interestVertices.size(); k++)
				{
					reducedGraphInterestVertices[k] = components[interestVertices[k]];
				}
				conditionalPMC(turnipInput);
				firstMomentSum += turnipInput.n * turnipInput.estimateFirstMoment;
			}
			firstMomentSum *= initialConditioningProb;
		}
		else
		{
			firstMomentSum = observations.size();
			if (variableMap.count("outputConditionalDistribution") > 0)
			{
				empiricalDistribution empiricalDist(false, nEdges, context);
				empiricalDist.hintDataCount(observations.size());
				for (std::vector<::networkReliability::subObs::basic>::iterator j = observations.begin(); j != observations.end(); j++)
				{
					empiricalDist.add(j->getState());
				}

				std::string outputConditionalDistribution = variableMap["outputConditionalDistribution"].as<std::string>();
				std::ofstream outputStream(outputConditionalDistribution.c_str(), std::ios_base::out);
				boost::archive::binary_oarchive outputArchive(outputStream, boost::archive::no_codecvt);
				outputArchive << empiricalDist;

				outputStream.flush();
				outputStream.close();
			}
		}
		mpfr_class estimate = firstMomentSum / totalSamples;
		std::cout << "Estimate is " << toString(estimate) << std::endl;
		return 0;
	}
}
int main(int argc, char **argv)
{
	return networkReliability::main(argc, argv);
}
