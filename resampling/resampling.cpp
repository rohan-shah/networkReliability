#include <boost/program_options.hpp>
#include "Arguments.h"
#include "NetworkReliabilityObs.h"
#include "NetworkReliabilitySubObs.h"
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
			("nPMC", boost::program_options::value<std::size_t>()->default_value(0ULL)->implicit_value(0ULL), "(int) Number of PMC samples to use")
			("outputConditionalDistribution", boost::program_options::value<std::string>(), "(path) File to output the empirical conditional distribution")
			("useSpatialDistances", boost::program_options::value<bool>()->default_value(false)->implicit_value(true), "(flag) Use spatial rather than combinatoric distances")
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
		context.setMinCut(true);

		std::size_t nPMC = variableMap["nPMC"].as<std::size_t>();

		int initialRadius;
		std::string message;
		if(!readInitialRadius(variableMap, initialRadius, message))
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
		std::vector<boost::accumulators::accumulator_set<calculation_type, boost::accumulators::stats<boost::accumulators::tag::sum> > > probabilities(initialRadius + 1, zeroInitialisedAccumulator);
		//mean conditioningProbability
		std::vector<boost::accumulators::accumulator_set<calculation_type, boost::accumulators::stats<boost::accumulators::tag::sum> > > conditioningProbabilities(initialRadius + 1, zeroInitialisedAccumulator);
		const std::vector<int>& interestVertices = context.getInterestVertices();

		//working data for algorithms
		std::vector<int> components;
		std::vector<boost::default_color_type> colorMap;
		boost::detail::depth_first_visit_restricted_impl_helper<Context::internalGraph>::stackType stack;

		calculation_type estimate = 0;

		//If the usePMC flag is set, don't use splitting on the last step. Instead use PMC. 
		int finalSplittingStep;
		if (nPMC > 0) finalSplittingStep = initialRadius - 1;
		else finalSplittingStep = initialRadius;

		//additional working data for getRadius1ReducedGraph
		std::vector<int> edgeCounts;
		std::vector<int> reducedGraphInterestVertices(interestVertices.size());
		ConditionalTurnipInput turnipInput(randomSource, NULL, reducedGraphInterestVertices);
		turnipInput.exponentialRate = -boost::multiprecision::log(mpfr_class(1 - opProbability));

		//If we specify the useMinCut option then we need to do resampling. This vector holds the probabilities
		std::vector<double> resamplingProbabilities;

		for (int i = 0; i < n; i++)
		{
			NetworkReliabilityObs currentObs = NetworkReliabilityObs::constructConditional(context, randomSource);
			conditioningProbabilities[0](currentObs.getConditioningProb());
			NetworkReliabilitySubObs subObs = currentObs.getSubObservation(initialRadius);

			if (isSingleComponent(context, subObs.getState(), components, stack, colorMap))
			{
				probabilities[0](0);
			}
			else
			{
				probabilities[0](1);
				observations.push_back(std::move(subObs));
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
				conditioningProbabilities[splittingLevel+1](newObs.getConditioningProb());
				NetworkReliabilitySubObs sub = newObs.getSubObservation(initialRadius - splittingLevel - 1);
				if(!isSingleComponent(context, sub.getState(), components, stack, colorMap))
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
			int previousSize = observations.size();
			//resampling step
			observations.clear();
			resamplingProbabilities.clear();
			mpfr_class sum = 0;
			for (std::vector<NetworkReliabilitySubObs>::iterator j = nextStepObservations.begin(); j != nextStepObservations.end(); j++)
			{
				sum += j->getConditioningProb();
				resamplingProbabilities.push_back(j->getConditioningProb().convert_to<double>());
			}
			if(sum.convert_to<double>() == 0) throw std::runtime_error("Sum of importance weights was zero");
			mpfr_class averageWeight = sum / previousSize;
			aliasMethod::aliasMethod alias(resamplingProbabilities, sum.convert_to<double>());
			for (int k = 0; k < previousSize; k++)
			{
				observations.push_back(nextStepObservations[alias(randomSource)].copyWithConditioningProb(averageWeight));
			}
		}
		if (nPMC > 0)
		{
			for (std::vector<NetworkReliabilitySubObs>::iterator j = observations.begin(); j != observations.end(); j++)
			{
				Context::internalGraph reducedGraph;
				j->getRadius1ReducedGraph(reducedGraph, turnipInput.minimumInoperative, edgeCounts, components, stack, colorMap);
				turnipInput.n = (int)nPMC;
				turnipInput.graph = &reducedGraph;
				const std::size_t nReducedVertices = boost::num_vertices(reducedGraph);
				const std::size_t nReducedEdges = boost::num_edges(reducedGraph);
				//If the graph is DEFINITELY disconnected, add a value of 1 and continue.
				if (nReducedEdges == 0)
				{
					bool alwaysConnected = true;
					for (int interestCounter = 1; interestCounter < interestVertices.size(); interestCounter++)
					{
						if (components[interestVertices[interestCounter]] != components[interestVertices[0]])
						{
							alwaysConnected = false;
							break;
						}
					}
					if (alwaysConnected) throw std::runtime_error("Internal error");
					estimate += j->getConditioningProb();
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
				for (int k = 0; k < interestVertices.size(); k++)
				{
					reducedGraphInterestVertices[k] = components[interestVertices[k]];
				}
				conditionalPMC(turnipInput);
				estimate += j->getGeneratedObservationConditioningProb() * turnipInput.estimateFirstMoment;
			}
			estimate /= n;
		}
		else
		{
			estimate = (boost::accumulators::sum(probabilities[finalSplittingStep]) / n);
			if (variableMap.count("outputConditionalDistribution") > 0)
			{
				const std::size_t nEdges = context.getNEdges();
				std::string outputConditionalDistribution = variableMap["outputConditionalDistribution"].as<std::string>();
				std::ofstream outputStream(outputConditionalDistribution.c_str(), std::ios_base::binary);
				outputStream.write((char*)&nEdges, sizeof(std::size_t));
				const std::size_t sampleSize = observations.size();
				outputStream.write((char*)&sampleSize, sizeof(std::size_t));
				//write by individual bits
				int nStoredBits = 0;
				unsigned int currentValue = 0;
				for (std::vector<NetworkReliabilitySubObs>::iterator j = observations.begin(); j != observations.end(); j++)
				{
					const EdgeState* state = j->getState();
					for(int k = 0; k < nEdges; k++)
					{
						nStoredBits++;
						currentValue<<=1;
						currentValue += ((state[k] & OP_MASK) > 0);
						if(nStoredBits == sizeof(int)*8) 
						{
							outputStream.write((char*)&currentValue, sizeof(int));
							currentValue = nStoredBits = 0;
						}
					}
				}
				if(nStoredBits != 0) outputStream.write((char*)&currentValue, sizeof(int));
				//Alternative is to write by whole bytes
				/*for (std::vector<NetworkReliabilitySubObs>::iterator j = observations.begin(); j != observations.end(); j++)
				{
					const EdgeState* state = j->getState();
					for(int k = 0; k < nEdges; k++)
					{
						outputStream << (char)((state[k] & OP_MASK) > 0);
					}
				}*/
				outputStream.flush();
				outputStream.close();
			}
		}
returnEstimate:
		std::cout << "Estimate is " << toString(estimate) << std::endl;
		return 0;
	}
}
int main(int argc, char **argv)
{
	return networkReliability::main(argc, argv);
}
