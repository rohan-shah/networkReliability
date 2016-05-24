#include <boost/program_options.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include "arguments.h"
#include "obs/withResampling.h"
#include "subObs/withResampling.h"
#include "networkReliabilityObsCollection.h"
#include <vector>
#include "graphAlgorithms.h"
#include <boost/random/bernoulli_distribution.hpp>
#include "argumentsMPFR.h"
#include "conditionalPMC.h"
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
#include "networkReliabilityObsTree.h"
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
			("nPMC", boost::program_options::value<std::size_t>()->default_value(0ULL)->implicit_value(0ULL), "(int) Number of PMC samples to use")
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

		context contextObj = context::emptyContext();
		if(!readContext(variableMap, contextObj, opProbability))
		{
			std::cout << "Unable to construct context object" << std::endl;
			return 0;
		}
		const std::vector<int>& interestVertices = contextObj.getInterestVertices();

		std::size_t nPMC = variableMap["nPMC"].as<std::size_t>();
		if(nPMC == 0)
		{
			std::cout << "Input nPMC cannot be zero" << std::endl;
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

		resamplingInput resamplingInputs(contextObj);
		resamplingInputs.shouldOutputTree = variableMap.count("outputTree") > 0;
		resamplingInputs.n = n;

		std::string message;
		if(!readThresholds(variableMap, resamplingInputs.thresholds, message))
		{
			std::cout << message << std::endl;
			return 0;
		}
		resamplingInputs.finalSplittingStep = (int)(resamplingInputs.thresholds.size() - 2);
		std::vector<::networkReliability::subObs::withResampling> observations;

		//working data for graph algorithms
		std::vector<int> components;
		std::vector<boost::default_color_type> colorMap;
		boost::detail::depth_first_visit_restricted_impl_helper<context::internalGraph>::stackType stack;

		mpfr_class estimate = 0;

		resamplingOutput resamplingOutputs(observations, randomSource, contextObj, resamplingInputs.thresholds);
		doResampling(resamplingInputs, resamplingOutputs);
		if(resamplingOutputs.zeroEstimate) goto returnEstimate;
		else
		{
			//Working data for call to getReducedGraph
			std::vector<int> edgeCounts;
			std::vector<int> reducedGraphInterestVertices(interestVertices.size());
			
			ConditionalTurnipInput turnipInput(randomSource, NULL, reducedGraphInterestVertices);
			turnipInput.exponentialRate = -boost::multiprecision::log(mpfr_class(1 - opProbability));
			for (std::vector<::networkReliability::subObs::withResampling>::iterator j = observations.begin(); j != observations.end(); j++)
			{
				context::internalGraph reducedGraph;
				turnipInput.minimumInoperative = j->getMinCut();
				j->getReducedGraph(reducedGraph, edgeCounts, components, stack, colorMap);
				turnipInput.n = (int)nPMC;
				turnipInput.graph = &reducedGraph;
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
					estimate += j->getConditioningProb();
					continue;
				}
				turnipInput.edges.clear();
				turnipInput.edges.resize(nReducedEdges);
				turnipInput.edgeCounts.resize(nReducedEdges);
				context::internalGraph::edge_iterator current, end;
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
				estimate += j->getGeneratedObservationConditioningProb() * turnipInput.estimateFirstMoment;
			}
			estimate /= n;
		}
returnEstimate:
		if(resamplingInputs.shouldOutputTree)
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
