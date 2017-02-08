#include "approximateZeroVarianceWORWithVarianceImpl.h"
#include <boost/random/random_number_generator.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include "conditionalPoissonSequential.h"
#include "sampford.h"
#include <boost/graph/undirected_dfs.hpp>
namespace networkReliability
{
	namespace approximateZeroVarianceWithVarianceImpl
	{
		struct approximateZeroVarianceScratch
		{
			std::vector<int> maxFlowResults;
			allPointsMaxFlow::allPointsMaxFlowScratch<context::internalDirectedGraph, int> allPointsMaxFlowScratch;
			std::vector<context::internalDirectedGraph::edge_descriptor> vertexPredecessor;
			std::vector<boost::default_color_type> colorVector;
			std::vector<int> distanceVector;
		};
		struct varianceGraphVertex
		{
			int indexWithinDesign;
			mutable ::sampling::mpfr_class accumulatedMean, accumulatedVariance;
			boost::shared_ptr<std::vector<::sampling::mpfr_class> > inclusionProbabilities;
			mpfr_class d;
		};
		typedef boost::adjacency_list<boost::listS, boost::vecS, boost::directedS, boost::property<boost::vertex_name_t, varianceGraphVertex> > varianceGraph;
		struct accumulationVisitor : public boost::default_dfs_visitor
		{
		public:
			void finish_vertex(const varianceGraph::vertex_descriptor& u, const varianceGraph& graph)
			{
				varianceGraph::out_edge_iterator current, end;
				boost::tie(current, end) = boost::out_edges(u, graph);
				if(current == end) return;

				const varianceGraphVertex& currentVertexInfo = boost::get(boost::vertex_name, graph, u);
				currentVertexInfo.accumulatedMean = currentVertexInfo.accumulatedVariance = 0.0;
				const varianceGraphVertex& firstChildVertex = boost::get(boost::vertex_name, graph, boost::target(*current, graph));
				boost::shared_ptr<std::vector<::sampling::mpfr_class> > inclusionProbabilities = firstChildVertex.inclusionProbabilities;
				for(; current != end; current++)
				{
					const varianceGraphVertex& childVertexInfo = boost::get(boost::vertex_name, graph, boost::target(*current, graph));
					varianceGraph::out_edge_iterator current2 = boost::out_edges(u, graph).first;
					currentVertexInfo.accumulatedMean += childVertexInfo.accumulatedMean / (*inclusionProbabilities)[childVertexInfo.indexWithinDesign];
					currentVertexInfo.accumulatedVariance += childVertexInfo.accumulatedVariance / (*inclusionProbabilities)[childVertexInfo.indexWithinDesign];
					for(; current2 != end; current2++)
					{
						const varianceGraphVertex& childVertexInfo2 = boost::get(boost::vertex_name, graph, boost::target(*current2, graph));
						if(childVertexInfo.indexWithinDesign == childVertexInfo2.indexWithinDesign)
						{
							currentVertexInfo.accumulatedVariance += ((1 - (*inclusionProbabilities)[childVertexInfo.indexWithinDesign]) / (*inclusionProbabilities)[childVertexInfo.indexWithinDesign]) * (childVertexInfo.accumulatedMean / (*inclusionProbabilities)[childVertexInfo.indexWithinDesign]) * (childVertexInfo2.accumulatedMean / (*inclusionProbabilities)[childVertexInfo2.indexWithinDesign]);
						}
						else
						{
							int N = (int)(*inclusionProbabilities).size();
							mpfr_class secondOrder;
							if(((*inclusionProbabilities)[childVertexInfo.indexWithinDesign] == 1 && (*inclusionProbabilities)[childVertexInfo2.indexWithinDesign] == 1) || childVertexInfo.d == 0)
							{
								secondOrder = 1;
							}
							else
							{
								secondOrder = (*inclusionProbabilities)[childVertexInfo.indexWithinDesign] * (*inclusionProbabilities)[childVertexInfo2.indexWithinDesign] - N * (*inclusionProbabilities)[childVertexInfo.indexWithinDesign] * (1 - (*inclusionProbabilities)[childVertexInfo.indexWithinDesign]) * (*inclusionProbabilities)[childVertexInfo2.indexWithinDesign] * (1 - (*inclusionProbabilities)[childVertexInfo2.indexWithinDesign])/((N - 1) * childVertexInfo.d);
							}
							currentVertexInfo.accumulatedVariance += (1 - (*inclusionProbabilities)[childVertexInfo.indexWithinDesign] * (*inclusionProbabilities)[childVertexInfo2.indexWithinDesign] / secondOrder) * (childVertexInfo.accumulatedMean / (*inclusionProbabilities)[childVertexInfo.indexWithinDesign]) * (childVertexInfo2.accumulatedMean / (*inclusionProbabilities)[childVertexInfo2.indexWithinDesign]);
						}
					}
				}
			}

		};
	}
	int getMinCut(std::vector<int>::iterator capacity, std::vector<int>::iterator residual, const context::internalDirectedGraph& graph, const context::internalGraph& undirectedGraph, const std::vector<int>& interestVertices, approximateZeroVarianceWithVarianceImpl::approximateZeroVarianceScratch& scratch)
	{
		std::size_t nVertices = boost::num_vertices(graph);
		scratch.vertexPredecessor.resize(nVertices);
		scratch.colorVector.resize(nVertices);
		scratch.distanceVector.resize(nVertices);

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
	void approximateZeroVarianceWORWithVariance(approximateZeroVarianceWORWithVarianceArgs& args)
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

		//The graph used to help estimate the variance. The initial vertex is the root. 
		approximateZeroVarianceWithVarianceImpl::varianceGraph varianceEstimationGraph;
		approximateZeroVarianceWithVarianceImpl::varianceGraph::vertex_descriptor rootVertex = boost::add_vertex(varianceEstimationGraph);
		std::vector<approximateZeroVarianceWithVarianceImpl::varianceGraph::vertex_descriptor> previousStepVertices, nextStepVertices;

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
		approximateZeroVarianceWithVarianceImpl::approximateZeroVarianceScratch scratch;
		//Get out the vector that holds the flow
		std::vector<int> residualCapacities(n * nEdges * 2), newResidualCapacities(n * nEdges * 2);
		std::vector<int> minCutSize(n), newMinCutSize(n);

		//Initialise with the two initial choices - The first edge can be up or down. 
		std::vector<mpfr_class> trueDensity, newTrueDensity, importanceDensity;
		std::vector<mpfr_class> newImportanceDensity, productInclusion, newProductInclusion;

		trueDensity.reserve(n);
		newTrueDensity.reserve(n);
		importanceDensity.reserve(n);
		newImportanceDensity.reserve(n);
		productInclusion.reserve(n);
		newProductInclusion.reserve(n);

		std::fill(states.begin(), states.begin()+2*nEdges, 1);
		states[0] = states[1] = 0;
		minCutSize[0] = getMinCut(states.begin(), residualCapacities.begin(), graph, undirectedGraph, interestVertices, scratch);

		std::fill(states.begin()+2*nEdges, states.begin()+4*nEdges, 1);
		states[2*nEdges] = states[2*nEdges + 1] = HIGH_CAPACITY;
		minCutSize[1] = getMinCut(states.begin()+2*nEdges, residualCapacities.begin()+2*nEdges, graph, undirectedGraph, interestVertices, scratch);
		if(minCutSize[1] < HIGH_CAPACITY)
		{
			//In this case there are two initial particles
			trueDensity.push_back(inopProbability);
			trueDensity.push_back(opProbability);
			mpfr_class minCutDownProb = cachedInopPowers[minCutSize[0]];
			mpfr_class minCutUpProb = cachedInopPowers[minCutSize[1]];
			mpfr_class qTilde = inopProbability * minCutDownProb;
			qTilde = qTilde / (qTilde  + opProbability * minCutUpProb);
			importanceDensity.push_back(qTilde);
			importanceDensity.push_back(1 - qTilde);

			productInclusion.push_back(1);
			productInclusion.push_back(1);

			//Now for the variance graph
			boost::shared_ptr<std::vector<::sampling::mpfr_class> > inclusion(new std::vector<::sampling::mpfr_class>(2, 1.0));

			approximateZeroVarianceWithVarianceImpl::varianceGraph::vertex_descriptor firstVertex = boost::add_vertex(varianceEstimationGraph);
			approximateZeroVarianceWithVarianceImpl::varianceGraphVertex& firstVertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, firstVertex);
			firstVertexInfo.inclusionProbabilities = inclusion;
			firstVertexInfo.indexWithinDesign = 0;
			boost::add_edge(rootVertex, firstVertex, varianceEstimationGraph);
			previousStepVertices.push_back(firstVertex);
			
			approximateZeroVarianceWithVarianceImpl::varianceGraph::vertex_descriptor secondVertex = boost::add_vertex(varianceEstimationGraph);
			approximateZeroVarianceWithVarianceImpl::varianceGraphVertex& secondVertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, secondVertex);
			secondVertexInfo.inclusionProbabilities = inclusion;
			secondVertexInfo.indexWithinDesign = 0;
			boost::add_edge(rootVertex, secondVertex, varianceEstimationGraph);
			previousStepVertices.push_back(secondVertex);
		}
		else
		{
			//In this case there is only one initial particle
			trueDensity.push_back(inopProbability);
			importanceDensity.push_back(1);
			productInclusion.push_back(1);

			//Add extra vertex to vertex estimation graph. 
			boost::shared_ptr<boost::numeric::ublas::matrix<::sampling::mpfr_class> > joint(new boost::numeric::ublas::matrix<::sampling::mpfr_class>(1, 1));
			std::fill((*joint).data().begin(), (*joint).data().end(), 1);
			boost::shared_ptr<std::vector<::sampling::mpfr_class> > inclusion(new std::vector<::sampling::mpfr_class>(1, 1.0));

			approximateZeroVarianceWithVarianceImpl::varianceGraph::vertex_descriptor singleVertexDescriptor = boost::add_vertex(varianceEstimationGraph);
			approximateZeroVarianceWithVarianceImpl::varianceGraphVertex& vertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, singleVertexDescriptor);
			vertexInfo.inclusionProbabilities = inclusion;
			vertexInfo.indexWithinDesign = 0;
			boost::add_edge(rootVertex, singleVertexDescriptor, varianceEstimationGraph);
			previousStepVertices.push_back(singleVertexDescriptor);
		}
		//The choices for sampling
		std::vector<choice> choices;
		args.estimate = 0;
		for(int edgeCounter = 1; edgeCounter < (int)nEdges; edgeCounter++)
		{
			//Construct choices
			choices.clear();
			newImportanceDensity.clear();
			for(int particleCounter = 0; particleCounter < (int)trueDensity.size(); particleCounter++)
			{
				if(minCutSize[particleCounter] == 0)
				{
					args.estimate += trueDensity[particleCounter] / productInclusion[particleCounter];
					approximateZeroVarianceWithVarianceImpl::varianceGraphVertex& vertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, previousStepVertices[particleCounter]);
					vertexInfo.accumulatedMean = trueDensity[particleCounter];
					vertexInfo.accumulatedVariance = 0;
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
			newTrueDensity.clear();
			newMinCutSize.clear();
			newProductInclusion.clear();
			importanceDensity.clear();
			nextStepVertices.clear();
			if(newImportanceDensity.size() <= n)
			{
				boost::shared_ptr<boost::numeric::ublas::matrix<::sampling::mpfr_class> > joint(new boost::numeric::ublas::matrix<::sampling::mpfr_class>(newImportanceDensity.size(), newImportanceDensity.size()));
				std::fill((*joint).data().begin(), (*joint).data().end(), 1);
				boost::shared_ptr<std::vector<::sampling::mpfr_class> > inclusion(new std::vector<::sampling::mpfr_class>(newImportanceDensity.size(), 1.0));
				for(int choiceCounter = 0; choiceCounter < (int)newImportanceDensity.size(); choiceCounter++)
				{
					int parentIndex = choices[choiceCounter].parentIndex;
					memcpy(&*(newStates.begin()+choiceCounter*2*nEdges), &*(states.begin()+parentIndex*2*nEdges), sizeof(int)*2*nEdges);
					memcpy(&*(newResidualCapacities.begin()+choiceCounter*2*nEdges), &*(residualCapacities.begin()+parentIndex*2*nEdges), sizeof(int)*2*nEdges);
					if(choices[choiceCounter].edgePresent)
					{
						newStates[choiceCounter*2*nEdges + 2*edgeCounter] = newStates[choiceCounter*2*nEdges + 2*edgeCounter + 1] = HIGH_CAPACITY;
						newTrueDensity.push_back(trueDensity[parentIndex]*opProbability);
					}
					else
					{
						newStates[choiceCounter*2*nEdges + 2*edgeCounter] = newStates[choiceCounter*2*nEdges + 2*edgeCounter + 1] = 0;
						newTrueDensity.push_back(trueDensity[parentIndex]*inopProbability);
					}
					newMinCutSize.push_back(choices[choiceCounter].minCutSize);
					newProductInclusion.push_back(productInclusion[parentIndex]);

					approximateZeroVarianceWithVarianceImpl::varianceGraph::vertex_descriptor newVertex = boost::add_vertex(varianceEstimationGraph);
					approximateZeroVarianceWithVarianceImpl::varianceGraphVertex& vertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, newVertex);
					vertexInfo.inclusionProbabilities = inclusion;
					vertexInfo.indexWithinDesign = choiceCounter;
					vertexInfo.accumulatedMean = vertexInfo.accumulatedVariance = 0.0;
					boost::add_edge(previousStepVertices[parentIndex], newVertex, varianceEstimationGraph);
					nextStepVertices.push_back(newVertex);
				}
				importanceDensity.swap(newImportanceDensity);
			}
			else
			{
				indices.clear();

				samplingArgs.weights.clear();
				for(int i = 0; i < (int)choices.size(); i++) samplingArgs.weights.push_back(newImportanceDensity[i] / productInclusion[choices[i].parentIndex]);
				
				sampling::sampfordFromParetoNaive(samplingArgs, args.randomSource);
				boost::shared_ptr<std::vector<::sampling::mpfr_class> > inclusionProbabilities(new std::vector<::sampling::mpfr_class>());
				inclusionProbabilities->swap(samplingArgs.rescaledWeights);

				mpfr_class d = 0;
				for(int i = 0; i < (int)samplingArgs.rescaledWeights.size(); i++) 
				{
					d += samplingArgs.rescaledWeights[i] * (1 - samplingArgs.rescaledWeights[i]);
				}
				
				int counter = 0;
				for(std::vector<int>::iterator i = indices.begin(); i != indices.end(); i++)
				{
					int parentIndex = choices[*i].parentIndex;
					memcpy(&*(newStates.begin()+counter*2*nEdges), &*(states.begin()+choices[*i].parentIndex*2*nEdges), sizeof(int)*2*nEdges);
					memcpy(&*(newResidualCapacities.begin()+counter*2*nEdges), &*(residualCapacities.begin()+parentIndex*2*nEdges), sizeof(int)*2*nEdges);
					if(choices[*i].edgePresent)
					{
						newStates[counter*2*nEdges + 2*edgeCounter] = newStates[counter*2*nEdges + 2*edgeCounter + 1] = HIGH_CAPACITY;
						newTrueDensity.push_back(trueDensity[parentIndex]*opProbability);
					}
					else
					{
						newStates[counter*2*nEdges + 2*edgeCounter] = newStates[counter*2*nEdges + 2*edgeCounter + 1] = 0;
						newTrueDensity.push_back(trueDensity[parentIndex]*inopProbability);
					}
					if(newTrueDensity.back() > 100) throw std::runtime_error("Internal error");
					newProductInclusion.push_back(productInclusion[parentIndex] * (*inclusionProbabilities)[*i]);
					
					approximateZeroVarianceWithVarianceImpl::varianceGraph::vertex_descriptor newVertex = boost::add_vertex(varianceEstimationGraph);
					approximateZeroVarianceWithVarianceImpl::varianceGraphVertex& vertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, newVertex);
					vertexInfo.indexWithinDesign = *i;
					vertexInfo.inclusionProbabilities = inclusionProbabilities;
					vertexInfo.d = d;
					vertexInfo.accumulatedMean = vertexInfo.accumulatedVariance = 0.0;
					boost::add_edge(previousStepVertices[parentIndex], newVertex, varianceEstimationGraph);
					nextStepVertices.push_back(newVertex);
					
					newMinCutSize.push_back(choices[*i].minCutSize);
					importanceDensity.push_back(newImportanceDensity[*i]);
					counter++;
				}
			}
			nextStepVertices.swap(previousStepVertices);
			productInclusion.swap(newProductInclusion);
			newTrueDensity.swap(trueDensity);
			newStates.swap(states);
			newMinCutSize.swap(minCutSize);
		}
		for(int i = 0; i < (int)trueDensity.size(); i++)
		{
			args.estimate += trueDensity[i] / productInclusion[i];
		}
		//Variance estimation
		for(int i = 0; i < (int)trueDensity.size(); i++)
		{
			approximateZeroVarianceWithVarianceImpl::varianceGraphVertex& vertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, previousStepVertices[i]);
			vertexInfo.accumulatedMean = trueDensity[i];
			vertexInfo.accumulatedVariance = 0;
		}
		typedef boost::color_traits<boost::default_color_type> Color;
		std::vector<boost::default_color_type> colourVector(boost::num_vertices(varianceEstimationGraph), Color::white());
		boost::iterator_property_map<std::vector<boost::default_color_type>::iterator, boost::identity_property_map> colourMap(colourVector.begin());
		boost::depth_first_visit(varianceEstimationGraph, rootVertex, approximateZeroVarianceWithVarianceImpl::accumulationVisitor(), colourMap);
		args.varianceEstimate = boost::get(boost::vertex_name, varianceEstimationGraph, rootVertex).accumulatedVariance;
	}
}
