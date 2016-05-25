#include "approximateZeroVarianceWORImpl.h"
#include <boost/random/random_number_generator.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/graph/stoer_wagner_min_cut.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
namespace networkReliability
{
	namespace approximateZeroVarianceImpl
	{
		struct approximateZeroVarianceScratch
		{
			std::vector<int> maxFlowResults;
			allPointsMaxFlow::allPointsMaxFlowScratch<context::internalDirectedGraph, int> allPointsMaxFlowScratch;
			std::vector<context::internalDirectedGraph::edge_descriptor> vertexPredecessor;
			std::vector<boost::default_color_type> colorVector;
			std::vector<int> distanceVector;
		};
	}
	int getMinCut(std::vector<int>::iterator capacity, std::vector<int>::iterator residual, const context::internalDirectedGraph& graph, const context::internalGraph& undirectedGraph, const std::vector<int>& interestVertices, approximateZeroVarianceImpl::approximateZeroVarianceScratch& scratch)
	{
		std::size_t nVertices = boost::num_vertices(graph);
		std::size_t nEdges = boost::num_edges(graph);
		scratch.vertexPredecessor.resize(nEdges);
		scratch.colorVector.resize(nEdges);
		scratch.distanceVector.resize(nEdges);

		//Are we looking at the all-terminal reliability problem?
		if(interestVertices.size() == nVertices)
		{
			contextImpl::constant_property_map_vertices_size_type<context::internalGraph::edge_descriptor, 1L> edgeWeights;

			//BOOST_AUTO(parities, boost::make_one_bit_color_map(num_vertices(*graph), get(boost::vertex_index, *graph)));
			return (int)boost::stoer_wagner_min_cut(undirectedGraph, edgeWeights);//, boost::parity_map(parities));
		}
		//or are we looking at the k-terminal reliability problem?
		else
		{
			const std::size_t nInterestVertices = interestVertices.size();
			//Use the all-points max flow
			if(interestVertices.size() * (interestVertices.size() - 1) / 2 > nVertices - 1)
			{
				scratch.maxFlowResults.resize(nVertices * nVertices, std::numeric_limits<int>::max());
				allPointsMaxFlow::allPointsMaxFlow<context::internalDirectedGraph, int>(scratch.maxFlowResults.begin(), capacity, graph, scratch.allPointsMaxFlowScratch);
				int minimum = std::numeric_limits<int>::max();
				for(std::size_t i = 0; i < nInterestVertices; i++)
				{
					for(std::size_t j = i+1; j < nInterestVertices; j++)
					{
						minimum = std::min(minimum, scratch.maxFlowResults[interestVertices[i] + interestVertices[j] * nVertices]);
					}
				}
				return minimum;
			}
			//Otherwise just call max-flow a bunch of times (the naive version)
			else
			{
				scratch.maxFlowResults.resize(nInterestVertices * nInterestVertices, std::numeric_limits<int>::max());
				typedef boost::property_map<context::internalDirectedGraph, boost::edge_index_t>::const_type edgeIndexMapType;
				typedef boost::property_map<context::internalDirectedGraph, boost::vertex_index_t>::const_type vertexIndexMapType;
				typedef boost::iterator_property_map<typename std::vector<int>::iterator, edgeIndexMapType> edgeCapacityMapType;
				typedef boost::iterator_property_map<typename std::vector<context::internalDirectedGraph::edge_descriptor>::iterator, vertexIndexMapType> vertexPredecessorMapType;
				typedef boost::iterator_property_map<typename std::vector<boost::default_color_type>::iterator, vertexIndexMapType> colorMapType;
				typedef boost::iterator_property_map<typename std::vector<int>::iterator, vertexIndexMapType> distanceMapType;

				edgeIndexMapType edgeIndexMap = boost::get(boost::edge_index, graph);
				vertexIndexMapType vertexIndexMap = boost::get(boost::vertex_index, graph);
				edgeCapacityMapType residualCapacityMap(residual, edgeIndexMap);
				edgeCapacityMapType edgeCapacityMap(capacity, edgeIndexMap);
				vertexPredecessorMapType vertexPredecessorMap(scratch.vertexPredecessor.begin(), vertexIndexMap);
				colorMapType colorMap(scratch.colorVector.begin(), vertexIndexMap);
				distanceMapType distanceMap(scratch.distanceVector.begin(), vertexIndexMap);

				for(std::size_t i = 0; i < nInterestVertices; i++)
				{
					for(std::size_t j = i+1; j < nInterestVertices; j++)
					{
						scratch.maxFlowResults[i + j * nInterestVertices] = scratch.maxFlowResults[j + i * nInterestVertices] = boost::boykov_kolmogorov_max_flow(graph, interestVertices[i], interestVertices[j], boost::residual_capacity_map(residualCapacityMap).capacity_map(edgeCapacityMap).predecessor_map(vertexPredecessorMap).color_map(colorMap).distance_map(distanceMap));
					}
				}
				return *std::min_element(scratch.maxFlowResults.begin(), scratch.maxFlowResults.end());
			}
		}
	}

	void approximateZeroVarianceWOR(approximateZeroVarianceWORArgs& args)
	{
		std::size_t n = args.n;
		if(n < 1)
		{
			throw std::runtime_error("Input n must be positive");
		}
		mpfr_class opProbability = args.contextObj.getOperationalProbability();
		mpfr_class inopProbability = 1 - opProbability;
		const std::size_t nEdges = args.contextObj.getNEdges();
		const std::vector<int>& interestVertices = args.contextObj.getInterestVertices();
		const context::internalDirectedGraph& graph = args.contextObj.getDirectedGraph();
		const context::internalGraph& undirectedGraph = args.contextObj.getGraph();
		
		if(interestVertices.size() > 2 && args.optimiseMinCut)
		{
			throw std::runtime_error("Can only specify option optimiseMinCut with 2-terminal reliability");
		}

		boost::random::uniform_01<float,float> uniformReal;
		//Vector used for mincut calculations
		std::vector<int> states(2*nEdges*n);
		//Sum over all the n simulations
		args.estimate = 0;

		//Cache powers of the inopProbability 
		boost::scoped_array<mpfr_class> cachedInopPowers(new mpfr_class[nEdges+1]);
		for(std::size_t i = 0; i < nEdges+1; i++)
		{
			cachedInopPowers[i] = boost::multiprecision::pow(inopProbability, i);
		}
		//Temporaries for calculating max flow values
		approximateZeroVarianceImpl::approximateZeroVarianceScratch scratch;
		//Get out the vector that holds the flow
		std::vector<int> residualCapacities(n * nEdges * 2);
		//Initialise with the two initial choices - The first edge can be up or down. 
		std::vector<mpfr_class> weights;
		weights.reserve(n);

		std::fill(states.begin(), states.begin()+2*nEdges, 1);
		states[0] = states[1] = 0;

		std::fill(states.begin()+2*nEdges, states.begin()+4*nEdges, 1);
		states[0] = states[1] = HIGH_CAPACITY;
		int initialMinCutSizeUp = getMinCut(states.begin(), residualCapacities.begin(), graph, undirectedGraph, interestVertices, scratch);
		if(initialMinCutSizeUp < HIGH_CAPACITY)
		{
			//In this case there are two initial particles
			weights.push_back(inopProbability);
			weights.push_back(opProbability);
		}
		else
		{
			//In this case there is only one initial particle
			weights.push_back(inopProbability);
		}
		//The choices for sampling
		std::vector<std::pair<int, bool> > choices;
		std::vector<bool> canReuseMinCut(n);
		std::vector<int> minCutSizeDown(n), minCutSizeUp(n);
		/*for(edgeCounter = 1; edgeCounter < nEdges; edgeCounter++)
		{
				bool canReuseMinCut = false;
				int minCutSizeDown, minCutSizeUp;
				for(edgeCounter = 1; edgeCounter < nEdges; edgeCounter++)
				{
					state[2*edgeCounter] = state[2*edgeCounter+1] = 0;
					if(!canReuseMinCut || !args.optimiseMinCut) minCutSizeDown = args.contextObj.getMinCut(state);
					minCutDownProb = cachedInopPowers[minCutSizeDown];

					state[2*edgeCounter] = state[2*edgeCounter+1] = HIGH_CAPACITY;
					if(!canReuseMinCut || !args.optimiseMinCut) minCutSizeUp = args.contextObj.getMinCut(state);
					if(minCutSizeUp >= HIGH_CAPACITY)
					{
						minCutUpProb = 0;
					}
					else minCutUpProb = cachedInopPowers[ minCutSizeUp];

					qTilde = inopProbability * minCutDownProb;
					qTilde = qTilde / (qTilde  + opProbability * minCutUpProb);
					float random = uniformReal(args.randomSource);
					if(random < qTilde)
					{
						state[2*edgeCounter] = state[2*edgeCounter+1] = 0; 
						currentLikelihoodRatio *= inopProbability / qTilde;
						if(minCutSizeDown == 0)
						{
							indicatorValue = 1;
							break;
						}
					}
					else 
					{
						state[2*edgeCounter] = state[2*edgeCounter+1] = HIGH_CAPACITY; 
						currentLikelihoodRatio *= opProbability / (1-qTilde);
						if(minCutSizeUp >= HIGH_CAPACITY)
						{
							indicatorValue = 0;
							break;
						}
					}
					canReuseMinCut = (residualCapacityVector[2*edgeCounter] == 1) && (residualCapacityVector[2*edgeCounter+1] == 1);
				}
			}
			args.estimateFirstMoment += indicatorValue * currentLikelihoodRatio;
			args.estimateSecondMoment += indicatorValue * currentLikelihoodRatio * currentLikelihoodRatio;
		}
		args.estimateFirstMoment /= args.n;
		args.estimateSecondMoment /= args.n;*/
	}
}
