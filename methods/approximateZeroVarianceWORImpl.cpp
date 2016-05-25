#include "approximateZeroVarianceWORImpl.h"
#include <boost/random/random_number_generator.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/graph/stoer_wagner_min_cut.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include "sampford.h"
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
		sampling::sampfordFromParetoNaiveArgs samplingArgs;
		//Temporaries for calculating max flow values
		approximateZeroVarianceImpl::approximateZeroVarianceScratch scratch;
		//Get out the vector that holds the flow
		std::vector<int> residualCapacities(n * nEdges * 2);
		std::vector<bool> canReuseMinCut(n);
		std::vector<int> minCutSize(n);

		//Initialise with the two initial choices - The first edge can be up or down. 
		std::vector<mpfr_class> weights;
		weights.reserve(n);

		std::fill(states.begin(), states.begin()+2*nEdges, 1);
		states[0] = states[1] = 0;
		minCutSize[0] = getMinCut(states.begin(), residualCapacities.begin(), graph, undirectedGraph, interestVertices, scratch);
		canReuseMinCut[0] = true;

		std::fill(states.begin()+2*nEdges, states.begin()+4*nEdges, 1);
		states[2*nEdges] = states[2*nEdges+1] = HIGH_CAPACITY;
		minCutSize[1] = getMinCut(states.begin()+2*nEdges, residualCapacities.begin()+2*nEdges, graph, undirectedGraph, interestVertices, scratch);
		if(minCutSize[1] < HIGH_CAPACITY)
		{
			//In this case there are two initial particles
			weights.push_back(inopProbability);
			weights.push_back(opProbability);
			canReuseMinCut[1] = (residualCapacities[2*nEdges] == 1) && (residualCapacities[2*nEdges+1] == 1);
		}
		else
		{
			//In this case there is only one initial particle
			weights.push_back(inopProbability);
		}
		//The choices for sampling
		std::vector<std::pair<int, bool> > choices;
		std::vector<mpfr_class> newWeights;
		for(int edgeCounter = 1; edgeCounter < (int)nEdges; edgeCounter++)
		{
			//Construct choices
			choices.clear();
			newWeights.clear();
			for(int particleCounter = 0; particleCounter < (int)weights.size(); particleCounter++)
			{
				states[2*nEdges*particleCounter + 2*edgeCounter] = states[2*nEdges*particleCounter + 2*edgeCounter+1] = 0;
				int minCutSizeDown;
				if(canReuseMinCut[particleCounter] && args.optimiseMinCut) minCutSizeDown = minCutSize[particleCounter];
				else minCutSizeDown = getMinCut(states.begin()+particleCounter*2*nEdges, residualCapacities.begin()+particleCounter*2*nEdges, graph, undirectedGraph, interestVertices, scratch);
				mpfr_class minCutDownProb = cachedInopPowers[minCutSizeDown];

				int minCutSizeUp;
				states[2*nEdges*particleCounter + 2*edgeCounter] = states[2*nEdges*particleCounter + 2*edgeCounter+1] = HIGH_CAPACITY;
				if(canReuseMinCut[particleCounter] || args.optimiseMinCut) minCutSizeUp = minCutSize[particleCounter];
				else
				{
					minCutSizeUp = getMinCut(states.begin()+particleCounter*2*nEdges, residualCapacities.begin()+particleCounter*2*nEdges, graph, undirectedGraph, interestVertices, scratch);
				}
				mpfr_class minCutUpProb;
				if(minCutSizeUp >= HIGH_CAPACITY)
				{
					minCutUpProb = 0;
					choices.push_back(std::make_pair(particleCounter, false));
					newWeights.push_back(weights[particleCounter]  / inopProbability);
				}
				else
				{
					minCutUpProb = cachedInopPowers[ minCutSizeUp];
					mpfr_class qTilde = inopProbability * minCutDownProb;
					qTilde = qTilde / (qTilde  + opProbability * minCutUpProb);
					choices.push_back(std::make_pair(particleCounter, false));
					choices.push_back(std::make_pair(particleCounter, true));
					newWeights.push_back(weights[particleCounter]  * qTilde / inopProbability);
					newWeights.push_back(weights[particleCounter]  * (1 - qTilde) / opProbability);
				}
			}
			/*args.estimateFirstMoment += indicatorValue * currentLikelihoodRatio;
			args.estimateSecondMoment += indicatorValue * currentLikelihoodRatio * currentLikelihoodRatio;*/
		}
		args.estimate /= args.n;
	}
}
