#include "approximateZeroVarianceWORImpl.h"
#include <boost/random/random_number_generator.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include "sampford.h"
namespace networkReliability
{
	namespace approximateZeroVarianceWORImpl
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
	int getMinCut(std::vector<int>::iterator capacity, std::vector<int>::iterator residual, const context::internalDirectedGraph& graph, const context::internalGraph& undirectedGraph, const std::vector<int>& interestVertices, approximateZeroVarianceWORImpl::approximateZeroVarianceScratch& scratch)
	{
		std::size_t nVertices = boost::num_vertices(graph);
		std::size_t nEdges = boost::num_edges(graph);
		scratch.vertexPredecessor.resize(nEdges);
		scratch.colorVector.resize(nEdges);
		scratch.distanceVector.resize(nEdges);

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
	void approximateZeroVarianceWOR(approximateZeroVarianceWORArgs& args)
	{
		struct choice
		{
			choice(int parentIndex, bool edgePresent, int minCutSize)
				:parentIndex(parentIndex), edgePresent(edgePresent), minCutSize(minCutSize)
			{}
			int parentIndex;
			bool edgePresent;
			int minCutSize;
		};
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
		
		boost::random::uniform_01<float,float> uniformReal;
		//Vector used for mincut calculations
		std::vector<int> states(2*nEdges*n), newStates(2*nEdges*n);

		//Cache powers of the inopProbability 
		boost::scoped_array<mpfr_class> cachedInopPowers(new mpfr_class[nEdges+1]);
		for(std::size_t i = 0; i < nEdges+1; i++)
		{
			cachedInopPowers[i] = boost::multiprecision::pow(inopProbability, i);
		}
		sampling::sampfordFromParetoNaiveArgs samplingArgs;
		samplingArgs.n = n;
		std::vector<int>& indices = samplingArgs.indices;
		//Temporaries for calculating max flow values
		approximateZeroVarianceWORImpl::approximateZeroVarianceScratch scratch;
		//Get out the vector that holds the flow
		std::vector<int> residualCapacities(n * nEdges * 2), newResidualCapacities(n * nEdges * 2);
		std::vector<int> minCutSize(n), newMinCutSize(n);

		//Initialise with the two initial choices - The first edge can be up or down. 
		std::vector<mpfr_class> weights, newWeights, importanceDensity;
		std::vector<mpfr_class>& newImportanceDensity = samplingArgs.weights;
		std::vector<mpfr_class>& rescaledWeights = samplingArgs.rescaledWeights;

		weights.reserve(n);
		newWeights.reserve(n);
		importanceDensity.reserve(n);

		std::fill(states.begin(), states.begin()+2*nEdges, 1);
		states[0] = states[1] = 0;
		minCutSize[0] = getMinCut(states.begin(), residualCapacities.begin(), graph, undirectedGraph, interestVertices, scratch);

		std::fill(states.begin()+2*nEdges, states.begin()+4*nEdges, 1);
		states[2*nEdges] = states[2*nEdges + 1] = HIGH_CAPACITY;
		minCutSize[1] = getMinCut(states.begin()+2*nEdges, residualCapacities.begin()+2*nEdges, graph, undirectedGraph, interestVertices, scratch);
		if(minCutSize[1] < HIGH_CAPACITY)
		{
			//In this case there are two initial particles
			weights.push_back(inopProbability);
			weights.push_back(opProbability);
			mpfr_class minCutDownProb = cachedInopPowers[minCutSize[0]];
			mpfr_class minCutUpProb = cachedInopPowers[minCutSize[1]];
			mpfr_class qTilde = inopProbability * minCutDownProb;
			qTilde = qTilde / (qTilde  + opProbability * minCutUpProb);
			importanceDensity.push_back(qTilde);
			importanceDensity.push_back(1 - qTilde);
		}
		else
		{
			//In this case there is only one initial particle
			weights.push_back(inopProbability);
			importanceDensity.push_back(1);
		}
		//The choices for sampling
		std::vector<choice> choices;
		args.estimate = 0;
		for(int edgeCounter = 1; edgeCounter < (int)nEdges; edgeCounter++)
		{
			//Construct choices
			choices.clear();
			newImportanceDensity.clear();
			for(int particleCounter = 0; particleCounter < (int)weights.size(); particleCounter++)
			{
				if(minCutSize[particleCounter] == 0)
				{
					args.estimate += weights[particleCounter];
					continue;
				}
				states[2*nEdges*particleCounter + 2*edgeCounter] = states[2*nEdges*particleCounter + 2*edgeCounter+1] = 0;
				int minCutSizeDown = getMinCut(states.begin()+particleCounter*2*nEdges, residualCapacities.begin()+particleCounter*2*nEdges, graph, undirectedGraph, interestVertices, scratch);
				mpfr_class minCutDownProb = cachedInopPowers[minCutSizeDown];

				states[2*nEdges*particleCounter + 2*edgeCounter] = states[2*nEdges*particleCounter + 2*edgeCounter+1] = HIGH_CAPACITY;
				int minCutSizeUp = getMinCut(states.begin()+particleCounter*2*nEdges, residualCapacities.begin()+particleCounter*2*nEdges, graph, undirectedGraph, interestVertices, scratch);
				mpfr_class minCutUpProb;
				if(minCutSizeUp >= HIGH_CAPACITY)
				{
					minCutUpProb = 0;
					choices.push_back(choice(particleCounter, false, minCutSizeDown));
					newImportanceDensity.push_back(importanceDensity[particleCounter]);
				}
				else
				{
					minCutUpProb = cachedInopPowers[minCutSizeUp];
					mpfr_class qTilde = inopProbability * minCutDownProb;
					qTilde = qTilde / (qTilde  + opProbability * minCutUpProb);
					choices.push_back(choice(particleCounter, false, minCutSizeDown));
					choices.push_back(choice(particleCounter, true, minCutSizeUp));
					newImportanceDensity.push_back(importanceDensity[particleCounter]  * qTilde);
					newImportanceDensity.push_back(importanceDensity[particleCounter]  * (1 - qTilde));
				}
			}
			newWeights.clear();
			newMinCutSize.clear();
			importanceDensity.clear();
			if(newImportanceDensity.size() <= n)
			{
				for(int choiceCounter = 0; choiceCounter < (int)newImportanceDensity.size(); choiceCounter++)
				{
					memcpy(&*(newStates.begin()+choiceCounter*2*nEdges), &*(states.begin()+choices[choiceCounter].parentIndex*2*nEdges), sizeof(int)*2*nEdges);
					memcpy(&*(newResidualCapacities.begin()+choiceCounter*2*nEdges), &*(residualCapacities.begin()+choices[choiceCounter].parentIndex*2*nEdges), sizeof(int)*2*nEdges);
					if(choices[choiceCounter].edgePresent)
					{
						newStates[choiceCounter*2*nEdges + 2*edgeCounter] = newStates[choiceCounter*2*nEdges + 2*edgeCounter + 1] = HIGH_CAPACITY;
						newWeights.push_back(weights[choices[choiceCounter].parentIndex]*opProbability);
					}
					else
					{
						newStates[choiceCounter*2*nEdges + 2*edgeCounter] = newStates[choiceCounter*2*nEdges + 2*edgeCounter + 1] = 0;
						newWeights.push_back(weights[choices[choiceCounter].parentIndex]*inopProbability);
					}
					newMinCutSize.push_back(choices[choiceCounter].minCutSize);
				}
				importanceDensity.swap(newImportanceDensity);
			}
			else
			{
				indices.clear();
				sampling::sampfordFromParetoNaive(samplingArgs, args.randomSource);
				int counter = 0;
				for(std::vector<int>::iterator i = indices.begin(); i != indices.end(); i++)
				{
					memcpy(&*(newStates.begin()+counter*2*nEdges), &*(states.begin()+choices[*i].parentIndex*2*nEdges), sizeof(int)*2*nEdges);
					memcpy(&*(newResidualCapacities.begin()+counter*2*nEdges), &*(residualCapacities.begin()+choices[*i].parentIndex*2*nEdges), sizeof(int)*2*nEdges);
					if(choices[*i].edgePresent)
					{
						newStates[counter*2*nEdges + 2*edgeCounter] = newStates[counter*2*nEdges + 2*edgeCounter + 1] = HIGH_CAPACITY;
						newWeights.push_back(weights[choices[*i].parentIndex]*opProbability / rescaledWeights[*i]);
					}
					else
					{
						newStates[counter*2*nEdges + 2*edgeCounter] = newStates[counter*2*nEdges + 2*edgeCounter + 1] = 0;
						newWeights.push_back(weights[choices[*i].parentIndex]*inopProbability / rescaledWeights[*i]);
					}
					if(newWeights.back() > 100) throw std::runtime_error("Internal error");
					newMinCutSize.push_back(choices[*i].minCutSize);
					//importanceDensity.push_back(1);
					importanceDensity.push_back(newImportanceDensity[*i] / rescaledWeights[*i]);
					counter++;
				}
			}
			newWeights.swap(weights);
			newStates.swap(states);
			newMinCutSize.swap(minCutSize);
		}
		for(std::vector<mpfr_class>::iterator i = weights.begin(); i != weights.end(); i++)
		{
			args.estimate += *i;
		}
	}
}
