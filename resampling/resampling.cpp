#include <boost/program_options.hpp>
#include <boost/math/special_functions/binomial.hpp>
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
#include "empiricalDistribution.h"
#include "depth_first_search_restricted.hpp"
#include "connected_components_restricted.hpp"
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/math/distributions.hpp>
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
			("nPMC", boost::program_options::value<std::size_t>()->default_value(0ULL)->implicit_value(0ULL), "(int) Number of PMC samples to use")
			("useCompleteEnumeration", boost::program_options::value<bool>()->default_value(false)->implicit_value(true), "(flag) Use complete enumeration in the last step")
			("outputConditionalDistribution", boost::program_options::value<std::string>(), "(path) File to output the empirical conditional distribution")
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
		const std::vector<int>& interestVertices = context.getInterestVertices();

		std::size_t nPMC = variableMap["nPMC"].as<std::size_t>();
		bool useCompleteEnumeration = variableMap["useCompleteEnumeration"].as<bool>();
		if(useCompleteEnumeration && nPMC > 0) 
		{
			std::cout << "Only one of nPMC and useCompleteEnumeration can be specified" << std::endl;
			return 0;
		}
		if(useCompleteEnumeration && interestVertices.size() > 2)
		{
			std::cout << "Complete enumeration only available for 2-terminal network reliability" << std::endl;
			return 0;
		}

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
		std::vector<int> resamplingIndices;

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

		//If the usePMC flag is set, don't use splitting on the last step. Instead use PMC. 
		int finalSplittingStep;
		if (nPMC > 0 || useCompleteEnumeration) finalSplittingStep = thresholds.size()-2;
		else finalSplittingStep = thresholds.size()-1;

		//If we specify the useMinCut option then we need to do resampling. This vector holds the probabilities
		std::vector<double> resamplingProbabilities;

		for (int i = 0; i < n; i++)
		{
			NetworkReliabilityObs currentObs = NetworkReliabilityObs::constructConditional(context, randomSource);
			NetworkReliabilitySubObs subObs = currentObs.getSubObservation(thresholds[0]);

			if (subObs.getMinCut() >= HIGH_CAPACITY)
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
			resamplingIndices.clear();
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
			//resampling step
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
			for (int k = 0; k < previousSize; k++)
			{
				int index = alias(randomSource);
				observations.push_back(nextStepObservations[index].copyWithGeneratedObservationConditioningProb(averageWeight));
				resamplingIndices.push_back(index);
			}
		}
		if (nPMC > 0)
		{
			//Working data for call to getRadius1ReducedGraph
			std::vector<int> edgeCounts;
			std::vector<int> reducedGraphInterestVertices(interestVertices.size());
			
			ConditionalTurnipInput turnipInput(randomSource, NULL, reducedGraphInterestVertices);
			turnipInput.exponentialRate = -boost::multiprecision::log(mpfr_class(1 - opProbability));
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
		else if(useCompleteEnumeration)
		{
			int tooManyEdgesCount = 0;
			estimate = 0;
			std::vector<EdgeState> edgeStates;
			NetworkReliabilitySubObs::getRadius1ReducedGraphNoSelfWithWeightsInput reducedGraphInput(interestVertices);
			//This is used to check the connected components of the reduced graph (not used in the call to getRadius1ReducedGraphNoSelfWithWeights
			std::vector<int> components2;
			//Similarly, this is used for the connected components of the reduced graph. 
			boost::detail::depth_first_visit_restricted_impl_helper<NetworkReliabilitySubObs::reducedGraphWithProbabilities>::stackType reducedGraphStack;
			//If we resample the same value multiple times above, only do the complete enumeration once. This isn't possible if we have only two thresholds (one non-zero, one zero) because there's no resampling done
			std::vector<mpfr_class> completeEnumerationValues(resamplingIndices.size(), -2);

			int counter = 0;
			int cachedValues = 0;
			for (std::vector<NetworkReliabilitySubObs>::iterator j = observations.begin(); j != observations.end(); j++, counter++)
			{
				if(resamplingIndices.size() > 0 && completeEnumerationValues[resamplingIndices[counter]] > -1)
				{
					mpfr_class currentEstimate = completeEnumerationValues[resamplingIndices[counter]];
					std::size_t nUnfixedEdges = 0;
					const EdgeState* state = j->getState();
					for(int i = 0; i < nEdges; i++)
					{
						if(state[i] & UNFIXED_MASK) nUnfixedEdges++;
					}
					boost::math::binomial_distribution<mpfr_class> targetDistribution(nUnfixedEdges, inopProbability);
					mpfr_class conditional = boost::math::cdf(boost::math::complement(targetDistribution, j->getMinCut()-1));
					estimate += j->getGeneratedObservationConditioningProb() * currentEstimate / conditional;
					cachedValues++;
					continue;
				}
				else if(j->getMinCut() > 0)
				{
					mpfr_class currentEstimate = 0;
					//get out the reduced graph
					j->getRadius1ReducedGraphNoSelfWithWeights(reducedGraphInput);
					std::size_t nReducedEdges = boost::num_edges(reducedGraphInput.outputGraph);
					std::size_t nReducedVertices = boost::num_vertices(reducedGraphInput.outputGraph);
					
					edgeStates.resize(nReducedEdges);
					components2.resize(nReducedVertices);
					colorMap.resize(nReducedVertices);
					int reducedVertex1 = reducedGraphInput.reducedInterestVertices[0];
					int reducedVertex2 = reducedGraphInput.reducedInterestVertices[1];
					EdgeState* edgeStatePtr = &(edgeStates[0]);
					//Exhaustive enumeration is only really going to be feasible if the number of edges in the reduced graph is small. 
					if(nReducedEdges < 20)
					{
						unsigned long maximumState = 1UL << nReducedEdges;
						for(unsigned long stateCounter = 0; stateCounter < maximumState; stateCounter++)
						{
							//expand out the bitmask. Initially specify that every reduced edge that is UP, is UP because all the component edges are up
							for(int edgeCounter = 0; edgeCounter < nReducedEdges; edgeCounter++)
							{
								if(stateCounter & (1LL << edgeCounter))
								{
									edgeStates[edgeCounter] = UNFIXED_OP;
								}
								else edgeStates[edgeCounter] = UNFIXED_INOP;
							}
							std::fill(colorMap.begin(), colorMap.end(), Color::white()); 
							boost::connected_components_restricted(reducedGraphInput.outputGraph, &(components2[0]), &(colorMap[0]), reducedGraphStack, &(edgeStates[0]));
							bool currentGraphConnected = components2[reducedVertex1] == components2[reducedVertex2];
							//if it's disconnected, we want to find the non-reduced edge configurations that give us that reduced edge configuration
							if(!currentGraphConnected)
							{
								mpfr_class currentPart = 1;
								NetworkReliabilitySubObs::reducedGraphWithProbabilities::edge_iterator current, end;
								boost::tie(current, end) = boost::edges(reducedGraphInput.outputGraph);
								for(int reducedEdgeCounter = 0; reducedEdgeCounter != nReducedEdges; reducedEdgeCounter++, current++)
								{
									if(edgeStates[reducedEdgeCounter] & UNFIXED_OP) currentPart *= boost::get(boost::edge_op_probability, reducedGraphInput.outputGraph, *current);
									else currentPart *= boost::get(boost::edge_inop_probability, reducedGraphInput.outputGraph, *current);
								}
								currentEstimate += currentPart;
							}
						}
						boost::math::binomial_distribution<mpfr_class> targetDistribution(reducedGraphInput.nUnreducedEdges, inopProbability);
						mpfr_class conditional = boost::math::cdf(boost::math::complement(targetDistribution, j->getMinCut()-1));
						estimate += j->getGeneratedObservationConditioningProb() * currentEstimate / conditional;
						if(resamplingIndices.size() > 1) completeEnumerationValues[resamplingIndices[counter]] = currentEstimate;
					}
					else
					{
						if(tooManyEdgesCount == 0)
						{
							std::ofstream outputStream("./tooManyEdges.dat");
							boost::archive::text_oarchive outputArchive(outputStream);
							writeNetworkReliabilitySubObs(outputArchive, *j);
						}
						tooManyEdgesCount++;
						NetworkReliabilityObs obs = j->getObservation(randomSource);
						if(!isSingleComponent(context, obs.getState(), components, stack, colorMap))
						{
							estimate += j->getGeneratedObservationConditioningProb();
						}
					}
				}
				else estimate += j->getGeneratedObservationConditioningProb();
			}
			estimate /= n;
			std::cout << tooManyEdgesCount << " graphs had too many edges for complete enumeration" << std::endl;
			std::cout << "Used cached values " << cachedValues << " times" << std::endl;
		}
		else
		{
			estimate = (boost::accumulators::sum(probabilities[finalSplittingStep]) / n);
			if (variableMap.count("outputConditionalDistribution") > 0)
			{
				empiricalDistribution outputDistributions(false, nEdges);
				for(std::vector<NetworkReliabilitySubObs>::iterator i = observations.begin(); i != observations.end(); i++)
				{
					outputDistributions.add(i->getState());
				}
				try
				{
					std::ofstream outputStream(variableMap["outputConditionalDistribution"].as<std::string>().c_str(), std::ios_base::binary);
					boost::archive::binary_oarchive archive(outputStream);
					archive << outputDistributions;
				}
				catch(std::runtime_error& err)
				{
					std::cout << "Error saving empirical distributions to file: " << err.what() << std::endl;
					return 0;
				}
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
