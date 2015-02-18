#include <boost/program_options.hpp>
#include "Arguments.h"
#include "subObs/withResampling.h"
#include "obs/withResampling.h"
#include "NetworkReliabilityObsCollection.h"
#include <vector>
#include "graphAlgorithms.h"
#include <boost/random/bernoulli_distribution.hpp>
#include "includeMPFR.h"
#include "ArgumentsMPFR.h"
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/special_functions.hpp>
#include <fstream>
#include <boost/random/uniform_real_distribution.hpp>
#include "empiricalDistribution.h"
#include "depth_first_search_restricted.hpp"
#include "connected_components_restricted.hpp"
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/math/distributions.hpp>
#include "NetworkReliabilityObsTree.h"
#include "commonOptions.h"
#include "resamplingCommon.h"
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
			("outputTree", boost::program_options::value<std::string>(), "(path) File to output simulation tree to")
			("useSpatialDistances", boost::program_options::value<std::vector<double> >()->multitoken(), "(float) Input spatial distances must consist of two or three numbers numbers; A maximum distance, an optional minimum distance and the number of steps to take.")
			("help", "Display this message");
		addCompensate(options);

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
			std::cout << "Complete enumeration only available for 2-terminal network reliability" << std::endl;
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

		bool compensateResampling = variableMap["compensateResampling"].as<bool>();
		resamplingInput resamplingInputs(context, compensateResampling);
		resamplingInputs.shouldOutputTree = variableMap.count("outputTree") > 0;
		resamplingInputs.n = n;

		std::string message;
		if(!readThresholds(variableMap, resamplingInputs.thresholds, message))
		{
			std::cout << message << std::endl;
			return 0;
		}
		resamplingInputs.finalSplittingStep = (int)(resamplingInputs.thresholds.size()-2);
		std::vector<::networkReliability::subObs::withResampling> observations;

		//working data for graph algorithms
		std::vector<int> components;
		std::vector<boost::default_color_type> colorMap;
		boost::detail::depth_first_visit_restricted_impl_helper<Context::internalGraph>::stackType stack;

		mpfr_class estimate = 0;

		resamplingOutput resamplingOutputs(observations, randomSource, context, resamplingInputs.thresholds);
		doResampling(resamplingInputs, resamplingOutputs);
		if(resamplingOutputs.zeroEstimate) goto returnEstimate;
		else
		{
			int tooManyEdgesCount = 0;
			std::vector<EdgeState> edgeStates;
			::networkReliability::subObs::subObs::getReducedGraphNoSelfWithWeightsInput reducedGraphInput(interestVertices);
			//This is used to check the connected components of the reduced graph (not used in the call to getReducedGraphNoSelfWithWeights
			std::vector<int> components2;
			//Similarly, this is used for the connected components of the reduced graph. 
			boost::detail::depth_first_visit_restricted_impl_helper<::networkReliability::subObs::subObs::reducedGraphWithProbabilities>::stackType reducedGraphStack;
			//Holds the subObservations that couldn't be completely enumerated
			NetworkReliabilityObsCollection tooManyEdges(&context, *(resamplingInputs.thresholds.rbegin()+1));

			int counter = 0;
			for (std::vector<::networkReliability::subObs::withResampling>::iterator j = observations.begin(); j != observations.end(); j++, counter++)
			{
				if(j->getMinCut() > 0)
				{
					mpfr_class currentEstimate = 0;
					//get out the reduced graph
					j->getReducedGraphNoSelfWithWeights(reducedGraphInput);
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
							for(std::size_t edgeCounter = 0; edgeCounter < nReducedEdges; edgeCounter++)
							{
								if(stateCounter & (1UL << edgeCounter))
								{
									edgeStatePtr[edgeCounter] = UNFIXED_OP;
								}
								else edgeStatePtr[edgeCounter] = UNFIXED_INOP;
							}
							std::fill(colorMap.begin(), colorMap.end(), Color::white()); 
							boost::connected_components_restricted(reducedGraphInput.outputGraph, &(components2[0]), &(colorMap[0]), reducedGraphStack, edgeStatePtr);
							bool currentGraphConnected = components2[reducedVertex1] == components2[reducedVertex2];
							//if it's disconnected, we want to find the non-reduced edge configurations that give us that reduced edge configuration
							if(!currentGraphConnected)
							{
								mpfr_class currentPart = 1;
								::networkReliability::subObs::subObs::reducedGraphWithProbabilities::edge_iterator current, end;
								boost::tie(current, end) = boost::edges(reducedGraphInput.outputGraph);
								for(std::size_t reducedEdgeCounter = 0; reducedEdgeCounter != nReducedEdges; reducedEdgeCounter++, current++)
								{
									if(edgeStatePtr[reducedEdgeCounter] & UNFIXED_OP) currentPart *= boost::get(boost::edge_op_probability, reducedGraphInput.outputGraph, *current);
									else currentPart *= boost::get(boost::edge_inop_probability, reducedGraphInput.outputGraph, *current);
								}
								currentEstimate += currentPart;
							}
						}
						boost::math::binomial_distribution<mpfr_class> targetDistribution(reducedGraphInput.nUnreducedEdges, inopProbability);
						mpfr_class conditional = boost::math::cdf(boost::math::complement(targetDistribution, j->getMinCut()-1));
						estimate += j->getGeneratedObservationConditioningProb() * currentEstimate / conditional;
					}
					else
					{
						tooManyEdges.add(*j);
						tooManyEdgesCount++;
						::networkReliability::obs::withResampling obs = ::networkReliability::subObs::getObservation<::networkReliability::subObs::withResampling>::get(*j, randomSource);
						if(!isSingleComponent(context, obs.getState(), components, stack, colorMap))
						{
							estimate += j->getGeneratedObservationConditioningProb();
						}
					}
				}
				else estimate += j->getGeneratedObservationConditioningProb();
			}
			if(tooManyEdgesCount > 0)
			{
				std::ofstream outputStream("./tooManyEdges.dat", std::ios_base::out | std::ios_base::binary);
				boost::archive::binary_oarchive outputArchive(outputStream, boost::archive::no_codecvt);
				outputArchive << tooManyEdges;
			}

			estimate /= n;
			std::cout << tooManyEdgesCount << " graphs had too many edges for complete enumeration" << std::endl;
		}
	returnEstimate:
		if(variableMap.count("outputTree") > 0)
		{
			std::cout << "Beginning tree layout....";
			bool success = resamplingOutputs.tree.layout();
			if(!success) std::cout << "Unable to lay out tree. Was graphviz support enabled? " << std::endl;
			else std::cout << "Done" << std::endl;
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
