#include "NetworkReliabilityObs.h"
#include <boost/program_options.hpp>
#include "Arguments.h"
#include "ArgumentsMPFR.h"
#include "graphAlgorithms.h"
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/sum_kahan.hpp>
namespace networkReliability
{
	int main(int argc, char **argv)
	{
		boost::program_options::options_description options("Usage");
		options.add_options()
			("gridGraph", boost::program_options::value<int>(), "(int) The dimension of the square grid graph to use. Incompatible with graphFile and torusGraph. ")
			("graphFile", boost::program_options::value<std::string>(), "(string) The path to a graphml file. Incompatible with gridGraph")
			("completeGraph", boost::program_options::value<int>(), "(int) The number of vertices of the complete graph to use. ")
			("opProbability", boost::program_options::value<double>(), "(float) The probability that an edge is operational. ")
			("n", boost::program_options::value<std::size_t>(), "(int) The number of simulations to perform. ")
			("seed", boost::program_options::value<int>(), "(int) The random seed used to generate the random graphs. ")
			("interestVertices", boost::program_options::value<std::vector<int> >()->multitoken(), "(int) The vertices of interest, that should be connected. ")
			("splitting", boost::program_options::value<int>(), "(int) Should we estimating splitting level probabilities?")
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
			std::cout << 
				"This program estimates the probability that the given graph is unreliable for the given vertices. That is, if edges are retained with a certain probability, what is the probability that the specified vertices are not all in the same connected component?\n\n"
			;
			std::cout << options << std::endl;
			return 0;
		}

		std::size_t n;
		if(!readN(variableMap, n))
		{
			return 0;
		}

		double probability;
		if(!readProbability(variableMap, probability))
		{
			return 0;
		}

		Context context = Context::emptyContext();
		if(!readContext(variableMap, context, probability))
		{
			return 0;
		}
		const std::size_t nEdges = context.getNEdges();

		boost::mt19937 randomSource;
		readSeed(variableMap, randomSource);

		boost::detail::depth_first_visit_restricted_impl_helper<Context::internalGraph>::stackType stack;

		std::vector<int> components;
		std::vector<int> interestComponents;
		const std::vector<int> interestVertices = context.getInterestVertices();

		std::vector<boost::default_color_type> colorMap;
		int splitting = 0;
		if(variableMap.count("splitting") && variableMap["splitting"].as<int>() > 0)
		{
			splitting = variableMap["splitting"].as<int>();
		}
		int countDisconnected = 0;
		std::vector<boost::accumulators::accumulator_set<int, boost::accumulators::stats<boost::accumulators::tag::sum, boost::accumulators::tag::count> > > conditioningProbabilities(splitting + 1);
		std::vector<boost::accumulators::accumulator_set<int, boost::accumulators::stats<boost::accumulators::tag::sum, boost::accumulators::tag::count> > > probabilities(splitting + 1);

		for(int i = 0; i < n; i++)
		{
			NetworkReliabilityObs obs(context, randomSource);
			if(!isSingleComponent(context, obs.getState(), components, stack, colorMap))
			{
				countDisconnected++;
			}			
			if(splitting > 0)
			{
				const EdgeState* obsState = obs.getState();
				int nDeactivated = 0;
				for(int k = 0; k < nEdges; k++)
				{
					if(obsState[k] & INOP_MASK) nDeactivated++;
				}
				//we need to iteratively update the bits that are fixed too. So we need a new state vector....
				boost::shared_array<EdgeState> updatedObsState(new EdgeState[nEdges]);
				memcpy(updatedObsState.get(), obsState, sizeof(EdgeState)*nEdges);

				std::size_t nextMinEdges = context.getMinCutEdges();
				for(int j = splitting; j >= 0; j--)
				{
					NetworkReliabilityObs obsWithFixed(context, updatedObsState, 0, 0);

					NetworkReliabilitySubObs subObs = obsWithFixed.getSubObservation(j);
					int subObsFixedInop = 0;
					const EdgeState* subObsState = subObs.getState();
					for(int k = 0; k < nEdges; k++)
					{
						if(subObsState[k] & (FIXED_INOP | NEW_FIXED_INOP)) subObsFixedInop++;
					}

					if(nDeactivated < nextMinEdges)
					{
						conditioningProbabilities[j](0);
						break;
					}
					else conditioningProbabilities[j](1);
					nextMinEdges = subObsFixedInop + subObs.getMinCut();
					if(isSingleComponent(context, subObs.getState(), components, stack, colorMap))
					{
						probabilities[j](0);
						break;
					}
					else
					{
						probabilities[j](1);
					}
					//anything that's fixed in subObs stays fixed
					for(int k = 0; k < nEdges; k++)
					{
						if(subObsState[k] & FIXED_MASK) updatedObsState[k] = subObsState[k];
					}
				}
			}
		}
		if(splitting > 0)
		{
			for(int i = splitting; i >= 0; i--)
			{
				std::cout << "Step " << (splitting - i) << ", conditioning on event " << boost::accumulators::sum(conditioningProbabilities[i]) << " / " <<  boost::accumulators::count(conditioningProbabilities[i]) << " = " << (float)boost::accumulators::sum(conditioningProbabilities[i]) / (float)boost::accumulators::count(conditioningProbabilities[i]) << std::endl;
				std::cout << "Step " << (splitting - i) << ", simulating event " << boost::accumulators::sum(probabilities[i]) << " / " <<  boost::accumulators::count(probabilities[i]) << " = " << (float)boost::accumulators::sum(probabilities[i]) / (float)boost::accumulators::count(probabilities[i]) << std::endl;
			}
		}
		std::cout << "Estimated unreliability probability was " << countDisconnected << " / " << n << " = " << (float)countDisconnected / (float)n << std::endl;;

		return 0;
	}
}
int main(int argc, char **argv)
{
	return networkReliability::main(argc, argv);
}