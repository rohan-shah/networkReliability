#include "approximateZeroVarianceWORWithVarianceImpl.h"
#include <boost/random/random_number_generator.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include "conditionalPoissonSequential.h"
#include "sampford.h"
#include <boost/graph/undirected_dfs.hpp>
#include <boost/random/uniform_01.hpp>
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
		struct particleData
		{
			std::vector<int> capacity, residual, downEdges;
		};
		struct varianceGraphVertex
		{
			int indexWithinDesign;
			mutable ::sampling::mpfr_class accumulatedMean, accumulatedVariance;
			boost::shared_ptr<std::vector<::sampling::mpfr_class> > inclusionProbabilities;
			boost::shared_ptr<boost::numeric::ublas::matrix<mpfr_class> > secondOrder;
			int edgeCounter;
			boost::shared_ptr<particleData> data;
		};
		typedef boost::adjacency_list<boost::listS, boost::vecS, boost::directedS, boost::property<boost::vertex_name_t, varianceGraphVertex> > varianceGraph;
		struct particle
		{
		public:
			boost::shared_ptr<particleData> parentData;
			boost::shared_ptr<particleData> ownedData;
			mpfr_class weight, importanceDensity, trueDensity;
			bool hasNextEdge;
			int minCutSize;
			int parentIndex;
			int edgeCounter;
			particle()
			{}
			particle(particle&& other)
				: parentData(std::move(other.parentData)), ownedData(std::move(other.ownedData)), weight(std::move(other.weight)), importanceDensity(std::move(other.importanceDensity)), trueDensity(other.trueDensity), hasNextEdge(other.hasNextEdge), minCutSize(other.minCutSize), parentIndex(other.parentIndex), edgeCounter(other.edgeCounter)
			{}
			particle& operator=(particle&& other)
			{
				parentData.swap(other.parentData);
				ownedData.swap(other.ownedData);
				weight = other.weight;
				importanceDensity = other.importanceDensity;
				trueDensity = other.trueDensity;
				hasNextEdge = other.hasNextEdge;
				minCutSize = other.minCutSize;
				parentIndex = other.parentIndex;
				edgeCounter = other.edgeCounter;
				return *this;
			}
		};
		double importanceSamplingEstimate(int edgeCounter,  boost::shared_ptr<particleData> data, const context& contextObj, boost::scoped_array<mpfr_class>& cachedInopPowers, boost::mt19937& randomSource)
		{
			mpfr_class opProbability = contextObj.getOperationalProbability();
			mpfr_class inopProbability = 1 - opProbability;
			const std::size_t nEdges = contextObj.getNEdges();
			mpfr_class currentLikelihoodRatio = 1;
			int indicatorValue = -1;
			boost::random::uniform_01<float,float> uniformReal;
			for(int edgeCounter2 = edgeCounter + 1; edgeCounter2 < (int)nEdges; edgeCounter2++)
			{
				data->capacity[2*edgeCounter2] = data->capacity[2*edgeCounter2+1] = 0;
				int minCutSizeDown = contextObj.getMinCut(data->capacity);
				mpfr_class minCutDownProb = cachedInopPowers[minCutSizeDown];

				data->capacity[2*edgeCounter2] = data->capacity[2*edgeCounter2+1] = HIGH_CAPACITY;
				int minCutSizeUp = contextObj.getMinCut(data->capacity);
				mpfr_class minCutUpProb;
				if(minCutSizeUp >= HIGH_CAPACITY)
				{
					minCutUpProb = 0;
				}
				else minCutUpProb = cachedInopPowers[minCutSizeUp];

				mpfr_class qTilde = inopProbability * minCutDownProb;
				qTilde = qTilde / (qTilde  + opProbability * minCutUpProb);
				float random = uniformReal(randomSource);
				if(random < qTilde)
				{
					data->capacity[2*edgeCounter2] = data->capacity[2*edgeCounter2+1] = 0; 
					currentLikelihoodRatio *= inopProbability / qTilde;
					if(minCutSizeDown == 0)
					{
						indicatorValue = 1;
						break;
					}
				}
				else 
				{
					data->capacity[2*edgeCounter2] = data->capacity[2*edgeCounter2+1] = HIGH_CAPACITY; 
					currentLikelihoodRatio *= opProbability / (1-qTilde);
					if(minCutSizeUp >= HIGH_CAPACITY)
					{
						indicatorValue = 0;
						break;
					}
				}
			}
			if(indicatorValue < 0)
			{
				 throw std::runtime_error("Internal error");
			}
			return indicatorValue*currentLikelihoodRatio.convert_to<double>();
		}
		struct accumulationVisitor : public boost::default_dfs_visitor
		{
		public:
			accumulationVisitor(const context& contextObj, boost::scoped_array<mpfr_class>& cachedInopPowers, boost::mt19937& randomSource)
				: contextObj(contextObj), cachedInopPowers(cachedInopPowers), randomSource(randomSource)
			{}
			accumulationVisitor(const accumulationVisitor& other)
				: contextObj(other.contextObj), cachedInopPowers(other.cachedInopPowers), randomSource(other.randomSource)
			{}
			const context& contextObj;
			boost::scoped_array<mpfr_class>& cachedInopPowers;
			boost::mt19937& randomSource;
			void finish_vertex(const varianceGraph::vertex_descriptor& u, const varianceGraph& graph)
			{
				const varianceGraphVertex& currentVertexInfo = boost::get(boost::vertex_name, graph, u);
				varianceGraph::out_edge_iterator current, end;
				boost::tie(current, end) = boost::out_edges(u, graph);
				if(current == end)
				{
/*					double firstBit = 1;
					for(int i = 0; i <= currentVertexInfo.edgeCounter; i++)
					{
						if(currentVertexInfo.data->capacity[2*i] == 0) firstBit *= inopProbability;
						else firstBit *= opProbability;
					}
					if(currentVertexInfo.accumulatedMean == 0)
					{
						double sum = 0;
						const int simulationCount = 2;
						for(int i = 0; i < simulationCount; i++)
						{
							sum += importanceSamplingEstimate(currentVertexInfo.edgeCounter, currentVertexInfo.data, contextObj, cachedInopPowers, randomSource);
						}
						currentVertexInfo.accumulatedMean = firstBit * sum / simulationCount;
					}
					else currentVertexInfo.accumulatedVariance = 0.0;*/
					return;	
				}
				else
				{
					/*varianceGraph::out_edge_iterator copied = current;
					copied++;
					if(copied == end)
					{
						return;
					}*/
				}

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
							currentVertexInfo.accumulatedVariance += ((1 - (*inclusionProbabilities)[childVertexInfo.indexWithinDesign])) * (childVertexInfo.accumulatedMean / (*inclusionProbabilities)[childVertexInfo.indexWithinDesign]) * (childVertexInfo2.accumulatedMean / (*inclusionProbabilities)[childVertexInfo2.indexWithinDesign]);
						}
						else
						{
							mpfr_class secondOrder = (*childVertexInfo.secondOrder.get())(childVertexInfo.indexWithinDesign, childVertexInfo2.indexWithinDesign);
							currentVertexInfo.accumulatedVariance += (1 - (*inclusionProbabilities)[childVertexInfo.indexWithinDesign] * (*inclusionProbabilities)[childVertexInfo2.indexWithinDesign] / secondOrder) * (childVertexInfo.accumulatedMean / (*inclusionProbabilities)[childVertexInfo.indexWithinDesign]) * (childVertexInfo2.accumulatedMean / (*inclusionProbabilities)[childVertexInfo2.indexWithinDesign]);
						}
					}
				}
			}

		};
	}
	int getMinCut(std::vector<int> capacity, std::vector<int> residual, const context::internalDirectedGraph& graph, const context::internalGraph& undirectedGraph, const std::vector<int>& interestVertices, approximateZeroVarianceWithVarianceImpl::approximateZeroVarianceScratch& scratch)
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
			allPointsMaxFlow::allPointsMaxFlow<context::internalDirectedGraph, int>(scratch.maxFlowResults.begin(), capacity.begin(), graph, scratch.allPointsMaxFlowScratch);
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
			edgeCapacityMapType residualCapacityMap(residual.begin(), edgeIndexMap);
			edgeCapacityMapType edgeCapacityMap(capacity.begin(), edgeIndexMap);
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
		sampling::conditionalPoissonSequentialArgs samplingArgs(true);
		samplingArgs.n = n;
		std::vector<int>& indices = samplingArgs.indices;
		//Temporaries for calculating max flow values
		approximateZeroVarianceWithVarianceImpl::approximateZeroVarianceScratch scratch;

		std::vector<mpfr_class>& newImportanceDensity = samplingArgs.weights;

		std::vector<approximateZeroVarianceWithVarianceImpl::particle> particles, newParticles;
		{
			approximateZeroVarianceWithVarianceImpl::particle missingEdgeParticle;
			missingEdgeParticle.ownedData.reset(new approximateZeroVarianceWithVarianceImpl::particleData());
			missingEdgeParticle.hasNextEdge = false;
			missingEdgeParticle.parentIndex = -1;
			missingEdgeParticle.ownedData->capacity.resize(2*nEdges, 1);
			missingEdgeParticle.ownedData->residual.resize(2*nEdges, 1);
			missingEdgeParticle.ownedData->capacity[0] = missingEdgeParticle.ownedData->capacity[1] = 0;
			missingEdgeParticle.ownedData->residual[0] = missingEdgeParticle.ownedData->residual[1] = 0;
			missingEdgeParticle.minCutSize = getMinCut(missingEdgeParticle.ownedData->capacity, missingEdgeParticle.ownedData->residual, graph, undirectedGraph, interestVertices, scratch);
			missingEdgeParticle.weight = inopProbability;
			missingEdgeParticle.trueDensity = inopProbability;
			missingEdgeParticle.ownedData->downEdges.push_back(0);

			approximateZeroVarianceWithVarianceImpl::particle withEdgeParticle;
			withEdgeParticle.ownedData.reset(new approximateZeroVarianceWithVarianceImpl::particleData());
			withEdgeParticle.hasNextEdge = true;
			withEdgeParticle.parentIndex = -1;
			withEdgeParticle.ownedData->capacity.resize(2*nEdges, 1);
			withEdgeParticle.ownedData->residual.resize(2*nEdges, 1);
			withEdgeParticle.ownedData->capacity[0] = withEdgeParticle.ownedData->capacity[1] = HIGH_CAPACITY;
			withEdgeParticle.ownedData->residual[0] = withEdgeParticle.ownedData->residual[1] = HIGH_CAPACITY;
			withEdgeParticle.minCutSize = getMinCut(withEdgeParticle.ownedData->capacity, withEdgeParticle.ownedData->residual, graph, undirectedGraph, interestVertices, scratch);
			withEdgeParticle.weight = opProbability;
			withEdgeParticle.trueDensity = opProbability;

			if(withEdgeParticle.minCutSize < HIGH_CAPACITY)
			{
				mpfr_class minCutDownProb = cachedInopPowers[missingEdgeParticle.minCutSize];
				mpfr_class minCutUpProb = cachedInopPowers[withEdgeParticle.minCutSize];
				mpfr_class qTilde = inopProbability * minCutDownProb;
				qTilde = qTilde / (qTilde  + opProbability * minCutUpProb);
	
				missingEdgeParticle.importanceDensity = qTilde;
				withEdgeParticle.importanceDensity = 1 - qTilde;
				//In this case there are two initial particles
				particles.emplace_back(std::move(missingEdgeParticle));
				particles.emplace_back(std::move(withEdgeParticle));

				boost::shared_ptr<std::vector<::sampling::mpfr_class> > inclusion(new std::vector<::sampling::mpfr_class>(2, 1.0));
				boost::shared_ptr<boost::numeric::ublas::matrix<mpfr_class> > secondOrder(new boost::numeric::ublas::matrix<mpfr_class>(2, 2, 1));
				approximateZeroVarianceWithVarianceImpl::varianceGraph::vertex_descriptor firstVertex = boost::add_vertex(varianceEstimationGraph);
				approximateZeroVarianceWithVarianceImpl::varianceGraphVertex& firstVertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, firstVertex);
				firstVertexInfo.edgeCounter = 0;
				firstVertexInfo.data = particles[0].ownedData;
				firstVertexInfo.inclusionProbabilities = inclusion;
				firstVertexInfo.secondOrder = secondOrder;
				firstVertexInfo.indexWithinDesign = 0;
				boost::add_edge(rootVertex, firstVertex, varianceEstimationGraph);
				previousStepVertices.push_back(firstVertex);
				
				approximateZeroVarianceWithVarianceImpl::varianceGraph::vertex_descriptor secondVertex = boost::add_vertex(varianceEstimationGraph);
				approximateZeroVarianceWithVarianceImpl::varianceGraphVertex& secondVertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, secondVertex);
				secondVertexInfo.edgeCounter = 0;
				secondVertexInfo.data = particles[1].ownedData;
				secondVertexInfo.inclusionProbabilities = inclusion;
				secondVertexInfo.secondOrder = secondOrder;
				secondVertexInfo.indexWithinDesign = 1;
				boost::add_edge(rootVertex, secondVertex, varianceEstimationGraph);
				previousStepVertices.push_back(secondVertex);
			}
			else
			{
				missingEdgeParticle.importanceDensity = 1;
				particles.emplace_back(std::move(missingEdgeParticle));

				boost::shared_ptr<std::vector<::sampling::mpfr_class> > inclusion(new std::vector<::sampling::mpfr_class>(1, 1.0));
				boost::shared_ptr<boost::numeric::ublas::matrix<mpfr_class> > secondOrder(new boost::numeric::ublas::matrix<mpfr_class>(1, 1, 1));
				approximateZeroVarianceWithVarianceImpl::varianceGraph::vertex_descriptor singleVertexDescriptor = boost::add_vertex(varianceEstimationGraph);
				approximateZeroVarianceWithVarianceImpl::varianceGraphVertex& vertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, singleVertexDescriptor);
				vertexInfo.edgeCounter = 0;
				vertexInfo.data = particles.back().ownedData;
				vertexInfo.inclusionProbabilities = inclusion;
				vertexInfo.secondOrder = secondOrder;
				vertexInfo.indexWithinDesign = 0;
				boost::add_edge(rootVertex, singleVertexDescriptor, varianceEstimationGraph);
				previousStepVertices.push_back(singleVertexDescriptor);
			}
		}
		args.estimate = 0;
		for(int edgeCounter = 1; edgeCounter < (int)nEdges; edgeCounter++)
		{
			newParticles.clear();
			for(int particleCounter = 0; particleCounter < (int)particles.size(); particleCounter++)
			{
				approximateZeroVarianceWithVarianceImpl::particle& currentParticle = particles[particleCounter];
				if(currentParticle.minCutSize == 0)
				{
					args.estimate += currentParticle.weight;
					approximateZeroVarianceWithVarianceImpl::varianceGraphVertex& vertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, previousStepVertices[particleCounter]);
					vertexInfo.accumulatedMean = currentParticle.trueDensity;
					vertexInfo.accumulatedVariance = 0;
					continue;
				}
				currentParticle.ownedData->capacity[2*edgeCounter] = currentParticle.ownedData->capacity[2*edgeCounter + 1] = 0;
				currentParticle.ownedData->residual[2*edgeCounter] = currentParticle.ownedData->residual[2*edgeCounter + 1] = 0;
				int minCutSizeDown = getMinCut(currentParticle.ownedData->capacity, currentParticle.ownedData->residual, graph, undirectedGraph, interestVertices, scratch);

				currentParticle.ownedData->capacity[2*edgeCounter] = currentParticle.ownedData->capacity[2*edgeCounter + 1] = HIGH_CAPACITY;
				currentParticle.ownedData->residual[2*edgeCounter] = currentParticle.ownedData->residual[2*edgeCounter + 1] = HIGH_CAPACITY;
				int minCutSizeUp = getMinCut(currentParticle.ownedData->capacity, currentParticle.ownedData->residual, graph, undirectedGraph, interestVertices, scratch);
				if(minCutSizeUp >= HIGH_CAPACITY)
				{
					approximateZeroVarianceWithVarianceImpl::particle newParticle;
					newParticle.parentData = particles[particleCounter].ownedData;
					newParticle.hasNextEdge = false;
					newParticle.weight = currentParticle.weight * inopProbability;
					newParticle.trueDensity = currentParticle.trueDensity * inopProbability;
					newParticle.importanceDensity = currentParticle.importanceDensity;
					newParticle.minCutSize = minCutSizeDown;
					newParticle.parentIndex = particleCounter;
					newParticles.emplace_back(std::move(newParticle));
				}
				else
				{
					mpfr_class minCutDownProb = cachedInopPowers[minCutSizeDown];
					mpfr_class minCutUpProb = cachedInopPowers[minCutSizeUp];
					mpfr_class qTilde = inopProbability * minCutDownProb;
					qTilde = qTilde / (qTilde  + opProbability * minCutUpProb);
					
					approximateZeroVarianceWithVarianceImpl::particle upParticle, downParticle;

					upParticle.parentData = currentParticle.ownedData;
					upParticle.hasNextEdge = true;
					upParticle.weight = currentParticle.weight * opProbability;
					upParticle.trueDensity = currentParticle.trueDensity * opProbability;
					upParticle.importanceDensity = currentParticle.importanceDensity * (1 - qTilde);
					upParticle.minCutSize = minCutSizeUp;

					downParticle.parentData = currentParticle.ownedData;
					downParticle.hasNextEdge = false;
					downParticle.weight = currentParticle.weight * inopProbability;
					downParticle.trueDensity = currentParticle.trueDensity * inopProbability;
					downParticle.importanceDensity = currentParticle.importanceDensity * qTilde;
					downParticle.minCutSize = minCutSizeDown;

					downParticle.parentIndex = upParticle.parentIndex = particleCounter;
					newParticles.emplace_back(std::move(downParticle));
					newParticles.emplace_back(std::move(upParticle));
				}
			}
			particles.clear();
			nextStepVertices.clear();
			if(newParticles.size() <= n)
			{
				boost::shared_ptr<std::vector<::sampling::mpfr_class> > inclusion(new std::vector<::sampling::mpfr_class>(newParticles.size(), 1.0));
				boost::shared_ptr<boost::numeric::ublas::matrix<mpfr_class> > secondOrder(new boost::numeric::ublas::matrix<mpfr_class>(newParticles.size(), newParticles.size(), 1));
				for(int i = 0; i < (int)newParticles.size(); i++)
				{
					particles.emplace_back(std::move(newParticles[i]));
					approximateZeroVarianceWithVarianceImpl::particle& movedParticle = particles.back();
					if(!(movedParticle.ownedData))
					{
						movedParticle.ownedData.reset(new approximateZeroVarianceWithVarianceImpl::particleData());
						movedParticle.ownedData->capacity =  movedParticle.parentData->capacity;
						movedParticle.ownedData->residual =  movedParticle.parentData->residual;
						movedParticle.ownedData->downEdges = movedParticle.parentData->downEdges;
						if(movedParticle.hasNextEdge)
						{
							movedParticle.ownedData->capacity[2*edgeCounter] = movedParticle.ownedData->capacity[2*edgeCounter+1] = movedParticle.ownedData->residual[2*edgeCounter] = movedParticle.ownedData->residual[2*edgeCounter+1] = HIGH_CAPACITY;
						}
						else
						{
							movedParticle.ownedData->capacity[2*edgeCounter] = movedParticle.ownedData->capacity[2*edgeCounter+1] = movedParticle.ownedData->residual[2*edgeCounter] = movedParticle.ownedData->residual[2*edgeCounter+1] = 0;
							movedParticle.ownedData->downEdges.push_back(edgeCounter);
						}
					}
					else if(!movedParticle.hasNextEdge)
					{
						movedParticle.ownedData->downEdges.push_back(edgeCounter);
					}

					approximateZeroVarianceWithVarianceImpl::varianceGraph::vertex_descriptor newVertex = boost::add_vertex(varianceEstimationGraph);
					approximateZeroVarianceWithVarianceImpl::varianceGraphVertex& vertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, newVertex);
					vertexInfo.edgeCounter = edgeCounter;
					vertexInfo.data = movedParticle.ownedData;
					vertexInfo.inclusionProbabilities = inclusion;
					vertexInfo.secondOrder = secondOrder;
					vertexInfo.indexWithinDesign = i;
					vertexInfo.accumulatedMean = vertexInfo.accumulatedVariance = 0.0;
					boost::add_edge(previousStepVertices[movedParticle.parentIndex], newVertex, varianceEstimationGraph);
					nextStepVertices.push_back(newVertex);
				}
			}
			else
			{
				newImportanceDensity.clear();
				for(std::vector<approximateZeroVarianceWithVarianceImpl::particle>::iterator i = newParticles.begin(); i != newParticles.end(); i++)
				{
					newImportanceDensity.push_back(i->importanceDensity);
				}
				indices.clear();
				sampling::conditionalPoissonSequential(samplingArgs, args.randomSource);
				boost::shared_ptr<boost::numeric::ublas::matrix<mpfr_class> > secondOrder(new boost::numeric::ublas::matrix<mpfr_class>(newParticles.size(), newParticles.size()));
				conditionalPoissonSecondOrderInclusionProbabilities(samplingArgs, samplingArgs.inclusionProbabilities, *secondOrder.get());

				boost::shared_ptr<std::vector<::sampling::mpfr_class> > inclusionProbabilities(new std::vector<::sampling::mpfr_class>());
				inclusionProbabilities->swap(samplingArgs.inclusionProbabilities);

				int counter = 0;
				for(std::vector<int>::iterator i = indices.begin(); i != indices.end(); i++)
				{
					particles.emplace_back(std::move(newParticles[*i]));
					approximateZeroVarianceWithVarianceImpl::particle& movedParticle = particles.back();
					if(!(movedParticle.ownedData))
					{
						movedParticle.ownedData.reset(new approximateZeroVarianceWithVarianceImpl::particleData());
						movedParticle.ownedData->capacity =  movedParticle.parentData->capacity;
						movedParticle.ownedData->residual =  movedParticle.parentData->residual;
						movedParticle.ownedData->downEdges = movedParticle.parentData->downEdges;
						if(movedParticle.hasNextEdge)
						{
							movedParticle.ownedData->capacity[2*edgeCounter] = movedParticle.ownedData->capacity[2*edgeCounter+1] = movedParticle.ownedData->residual[2*edgeCounter] = movedParticle.ownedData->residual[2*edgeCounter+1] = HIGH_CAPACITY;
						}
						else
						{
							movedParticle.ownedData->capacity[2*edgeCounter] = movedParticle.ownedData->capacity[2*edgeCounter+1] = movedParticle.ownedData->residual[2*edgeCounter] = movedParticle.ownedData->residual[2*edgeCounter+1] = 0;
							movedParticle.ownedData->downEdges.push_back(edgeCounter);
						}
					}
					else if(!movedParticle.hasNextEdge)
					{
						movedParticle.ownedData->downEdges.push_back(edgeCounter);
					}
					movedParticle.importanceDensity /= (*inclusionProbabilities)[*i];
					movedParticle.weight /= (*inclusionProbabilities)[*i];
					
					approximateZeroVarianceWithVarianceImpl::varianceGraph::vertex_descriptor newVertex = boost::add_vertex(varianceEstimationGraph);
					approximateZeroVarianceWithVarianceImpl::varianceGraphVertex& vertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, newVertex);
					vertexInfo.edgeCounter = edgeCounter;
					vertexInfo.data = movedParticle.ownedData;
					vertexInfo.indexWithinDesign = *i;
					vertexInfo.inclusionProbabilities = inclusionProbabilities;
					vertexInfo.secondOrder = secondOrder;
					vertexInfo.accumulatedMean = vertexInfo.accumulatedVariance = 0.0;
					boost::add_edge(previousStepVertices[movedParticle.parentIndex], newVertex, varianceEstimationGraph);
					nextStepVertices.push_back(newVertex);
					
					counter++;
				}
			}
			nextStepVertices.swap(previousStepVertices);
		}
		for(int i = 0; i < (int)particles.size(); i++)
		{
			args.estimate += particles[i].weight;
		}
		//Variance estimation
		for(int i = 0; i < (int)particles.size(); i++)
		{
			approximateZeroVarianceWithVarianceImpl::varianceGraphVertex& vertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, previousStepVertices[i]);
			vertexInfo.accumulatedMean = particles[i].trueDensity;
			vertexInfo.accumulatedVariance = 0;
		}
		typedef boost::color_traits<boost::default_color_type> Color;
		std::vector<boost::default_color_type> colourVector(boost::num_vertices(varianceEstimationGraph), Color::white());
		boost::iterator_property_map<std::vector<boost::default_color_type>::iterator, boost::identity_property_map> colourMap(colourVector.begin());
		boost::depth_first_visit(varianceEstimationGraph, rootVertex, approximateZeroVarianceWithVarianceImpl::accumulationVisitor(args.contextObj, cachedInopPowers, args.randomSource), colourMap);
		args.varianceEstimate = boost::get(boost::vertex_name, varianceEstimationGraph, rootVertex).accumulatedVariance;
	}
}
