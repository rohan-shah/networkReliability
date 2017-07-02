#include "approximateZeroVarianceWORMergeWithVarianceImpl.h"
#include <boost/random/random_number_generator.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include "conditionalPoissonSequential.h"
#include "connected_components_fixed.hpp"
namespace networkReliability
{
	namespace approximateZeroVarianceWORMergeWithVarianceImpl
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
		public:
			varianceGraphVertex()
				: indexWithinDesign(-1), samplingStage(-1), mergedCount(-1), mergedProduct(-1), accumulatedMean(0)
			{}
			int indexWithinDesign;
			int samplingStage;
			int mergedCount;
			int mergedProduct;
			mutable ::sampling::mpfr_class accumulatedMean;
		};
		typedef boost::adjacency_list<boost::listS, boost::vecS, boost::directedS, boost::property<boost::vertex_name_t, varianceGraphVertex>, boost::property<boost::edge_name_t, bool> > varianceGraph;
		struct particle
		{
		public:
			boost::shared_ptr<particleData> parentData;
			boost::shared_ptr<particleData> ownedData;
			mpfr_class weight, importanceDensity;
			bool hasNextEdge;
			int minCutSize;
			int mergedCount;
			int mergedProduct;
			std::vector<int> parentIndices;
			std::vector<bool> hasNextEdgeVector;
			particle()
			{}
			particle(particle&& other)
				: parentData(std::move(other.parentData)), ownedData(std::move(other.ownedData)), weight(std::move(other.weight)), importanceDensity(std::move(other.importanceDensity)), hasNextEdge(other.hasNextEdge), minCutSize(other.minCutSize), mergedCount(other.mergedCount), mergedProduct(other.mergedProduct), parentIndices(std::move(other.parentIndices)), hasNextEdgeVector(std::move(other.hasNextEdgeVector))
			{}
			particle& operator=(particle&& other)
			{
				parentData.swap(other.parentData);
				ownedData.swap(other.ownedData);
				weight = other.weight;
				importanceDensity = other.importanceDensity;
				hasNextEdge = other.hasNextEdge;
				minCutSize = other.minCutSize;
				mergedCount = other.mergedCount;
				mergedProduct = mergedProduct;
				parentIndices.swap(other.parentIndices);
				hasNextEdgeVector.swap(other.hasNextEdgeVector);
				return *this;
			}
			bool order(const particle& other, int edgeCounter) const 
			{
				if(ownedData && other.ownedData)
				{
					return memcmp(&(ownedData->capacity[0]), &(other.ownedData->capacity[0]), sizeof(int)*2*edgeCounter) < 0;
				}
				else if(ownedData && !other.ownedData)
				{
					int comp = memcmp(&(ownedData->capacity[0]), &(other.parentData->capacity[0]), sizeof(int)*2*(edgeCounter - 1));
					if(comp != 0) return comp < 0;
					if(other.hasNextEdge == hasNextEdge) return false;
					if(other.hasNextEdge) return true;
					return false;
				}
				else if(!ownedData && other.ownedData)
				{
					int comp = memcmp(&(parentData->capacity[0]), &(other.ownedData->capacity[0]), sizeof(int)*2*(edgeCounter - 1));
					if(comp != 0) return comp < 0;
					if(other.hasNextEdge == hasNextEdge) return false;
					if(other.hasNextEdge) return true;
					return false;
				}
				else
				{
					int comp = memcmp(&(parentData->capacity[0]), &(other.parentData->capacity[0]), sizeof(int)*2*(edgeCounter - 1));
					if(comp != 0) return comp < 0;
					if(other.hasNextEdge == hasNextEdge) return false;
					if(other.hasNextEdge) return true;
					return false;
				}
			}
			bool matches(const particle& other, int edgeCounter) const 
			{
				if(ownedData && other.ownedData)
				{
					return memcmp(&(ownedData->capacity[0]), &(other.ownedData->capacity[0]), sizeof(int)*2*(edgeCounter+1)) == 0;
				}
				else if(ownedData && !other.ownedData)
				{
					if(memcmp(&(ownedData->capacity[0]), &(other.parentData->capacity[0]), sizeof(int)*2*edgeCounter) != 0) return false;
					return other.hasNextEdge == hasNextEdge;
				}
				else if(!ownedData && other.ownedData)
				{
					if(memcmp(&(other.ownedData->capacity[0]), &(parentData->capacity[0]), sizeof(int)*2*edgeCounter) != 0) return false;
					return other.hasNextEdge == hasNextEdge;
				}
				else
				{
					return hasNextEdge == other.hasNextEdge && memcmp(&(parentData->capacity[0]), &(other.parentData->capacity[0]), sizeof(int)*edgeCounter*2) == 0;
				}
			}
		};
	}
	int getMinCut(std::vector<int>& capacity, std::vector<int>& residual, const context::internalDirectedGraph& graph, const context::internalGraph& undirectedGraph, const std::vector<int>& interestVertices, approximateZeroVarianceWORMergeWithVarianceImpl::approximateZeroVarianceScratch& scratch)
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

			std::fill(residual.begin(), residual.end(), 0);
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
	void approximateZeroVarianceWORMergeWithVariance(approximateZeroVarianceWORMergeWithVarianceArgs& args)
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
		const std::size_t nVertices = boost::num_vertices(graph);

		boost::detail::depth_first_visit_fixed_impl_helper<context::internalGraph>::stackType componentsStack;
		std::vector<int> components(nVertices);
 
		std::vector<std::pair<int, int> > edges(nEdges);
		{
			context::internalGraph::edge_iterator current, end;
			boost::tie(current, end) = boost::edges(undirectedGraph);
			for(; current != end; current++)
			{
				int index = boost::get(boost::edge_index, undirectedGraph, *current), target = boost::target(*current, undirectedGraph), source = boost::source(*current, undirectedGraph);
				edges[index] = std::make_pair(target, source);
				
			}
		}
		
		boost::random::uniform_01<float,float> uniformReal;

		///The graph used to help estimate the variance. The initial vertex is the root.
		approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraph varianceEstimationGraph;
		boost::numeric::ublas::matrix<int> graphVertices(nEdges, n, -1);
		approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraph::vertex_descriptor rootVertex = boost::add_vertex(varianceEstimationGraph);

		std::vector<std::vector<::sampling::mpfr_class> > allInclusionProbabilities;
		std::vector<boost::numeric::ublas::matrix<mpfr_class> > allSecondOrderInclusionProbabilities;

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
		approximateZeroVarianceWORMergeWithVarianceImpl::approximateZeroVarianceScratch scratch;
		scratch.colorVector.resize(nVertices);

		//Initialise with the two initial choices - The first edge can be up or down. 
		std::vector<approximateZeroVarianceWORMergeWithVarianceImpl::particle> particles, newParticles;
		{
			approximateZeroVarianceWORMergeWithVarianceImpl::particle missingEdgeParticle;
			missingEdgeParticle.ownedData.reset(new approximateZeroVarianceWORMergeWithVarianceImpl::particleData());
			missingEdgeParticle.hasNextEdge = false;
			missingEdgeParticle.ownedData->capacity.resize(2*nEdges, 1);
			missingEdgeParticle.ownedData->residual.resize(2*nEdges, 1);
			missingEdgeParticle.ownedData->capacity[0] = missingEdgeParticle.ownedData->capacity[1] = 0;
			missingEdgeParticle.ownedData->residual[0] = missingEdgeParticle.ownedData->residual[1] = 0;
			missingEdgeParticle.minCutSize = getMinCut(missingEdgeParticle.ownedData->capacity, missingEdgeParticle.ownedData->residual, graph, undirectedGraph, interestVertices, scratch);
			missingEdgeParticle.weight = inopProbability;
			missingEdgeParticle.ownedData->downEdges.push_back(0);
			missingEdgeParticle.mergedCount = missingEdgeParticle.mergedProduct = 1;

			approximateZeroVarianceWORMergeWithVarianceImpl::particle withEdgeParticle;
			withEdgeParticle.ownedData.reset(new approximateZeroVarianceWORMergeWithVarianceImpl::particleData());
			withEdgeParticle.hasNextEdge = true;
			withEdgeParticle.ownedData->capacity.resize(2*nEdges, 1);
			withEdgeParticle.ownedData->residual.resize(2*nEdges, 1);
			withEdgeParticle.ownedData->capacity[0] = withEdgeParticle.ownedData->capacity[1] = HIGH_CAPACITY;
			withEdgeParticle.ownedData->residual[0] = withEdgeParticle.ownedData->residual[1] = HIGH_CAPACITY;
			withEdgeParticle.minCutSize = getMinCut(withEdgeParticle.ownedData->capacity, withEdgeParticle.ownedData->residual, graph, undirectedGraph, interestVertices, scratch);
			withEdgeParticle.weight = opProbability;
			withEdgeParticle.mergedCount = withEdgeParticle.mergedProduct = 1;

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

				//Code relating to the graph.
				std::vector<::sampling::mpfr_class> inclusion(2, 1.0);
				allInclusionProbabilities.emplace_back(std::move(inclusion));
				boost::numeric::ublas::matrix<mpfr_class> secondOrder(2, 2, 1);
				allSecondOrderInclusionProbabilities.emplace_back(std::move(secondOrder));

				approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraph::vertex_descriptor firstVertex = boost::add_vertex(varianceEstimationGraph);
				approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraphVertex& firstVertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, firstVertex);
				firstVertexInfo.samplingStage = 0;
				firstVertexInfo.indexWithinDesign = 0;
				firstVertexInfo.mergedCount = firstVertexInfo.mergedProduct = 1;
				boost::add_edge(rootVertex, firstVertex, varianceEstimationGraph);
				graphVertices(0, 0) = (int)firstVertex;

				approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraph::vertex_descriptor secondVertex = boost::add_vertex(varianceEstimationGraph);
				approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraphVertex& secondVertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, secondVertex);
				secondVertexInfo.samplingStage = 0;
				secondVertexInfo.indexWithinDesign = 1;
				secondVertexInfo.mergedCount = secondVertexInfo.mergedProduct = 1;
				boost::add_edge(rootVertex, secondVertex, varianceEstimationGraph);
				graphVertices(0, 1) = (int)secondVertex;
			}
			else
			{
				missingEdgeParticle.importanceDensity = 1;
				particles.emplace_back(std::move(missingEdgeParticle));

				std::vector<::sampling::mpfr_class> inclusion(1, 1.0);
				allInclusionProbabilities.emplace_back(std::move(inclusion));
				boost::numeric::ublas::matrix<mpfr_class> secondOrder(1, 1, 1);
				allSecondOrderInclusionProbabilities.emplace_back(std::move(secondOrder));

				approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraph::vertex_descriptor vertex = boost::add_vertex(varianceEstimationGraph);
				approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraphVertex& vertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, vertex);
				vertexInfo.samplingStage = 0;
				vertexInfo.indexWithinDesign = 0;
				vertexInfo.mergedCount = vertexInfo.mergedProduct = 1;
				boost::add_edge(rootVertex, vertex, varianceEstimationGraph);
				graphVertices(0, 0) = (int)vertex;
			}
		}
		args.estimate = 0;
		for(int edgeCounter = 1; edgeCounter < (int)nEdges; edgeCounter++)
		{
			newParticles.clear();
			for(int particleCounter = 0; particleCounter < (int)particles.size(); particleCounter++)
			{
				approximateZeroVarianceWORMergeWithVarianceImpl::particle& currentParticle = particles[particleCounter];
				if(currentParticle.minCutSize == 0)
				{
					args.estimate += currentParticle.weight;
					approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraphVertex& vertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, graphVertices(edgeCounter-1, particleCounter));
					vertexInfo.accumulatedMean = 1;
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
					approximateZeroVarianceWORMergeWithVarianceImpl::particle newParticle;
					newParticle.parentData = particles[particleCounter].ownedData;
					newParticle.hasNextEdge = false;
					newParticle.hasNextEdgeVector.push_back(false);
					newParticle.weight = currentParticle.weight * inopProbability;
					newParticle.importanceDensity = currentParticle.importanceDensity;
					newParticle.minCutSize = minCutSizeDown;
					newParticle.parentIndices.push_back(particleCounter);
					newParticles.emplace_back(std::move(newParticle));
					newParticle.mergedCount = currentParticle.mergedCount;
					newParticle.mergedProduct = 1;
				}
				else
				{
					mpfr_class minCutDownProb = cachedInopPowers[minCutSizeDown];
					mpfr_class minCutUpProb = cachedInopPowers[minCutSizeUp];
					mpfr_class qTilde = inopProbability * minCutDownProb;
					qTilde = qTilde / (qTilde  + opProbability * minCutUpProb);
					
					approximateZeroVarianceWORMergeWithVarianceImpl::particle upParticle, downParticle;

					upParticle.parentData = currentParticle.ownedData;
					upParticle.hasNextEdge = true;
					upParticle.hasNextEdgeVector.push_back(true);
					upParticle.parentIndices.push_back(particleCounter);
					upParticle.weight = currentParticle.weight * opProbability;
					upParticle.importanceDensity = currentParticle.importanceDensity * (1 - qTilde);
					upParticle.minCutSize = minCutSizeUp;
					upParticle.mergedCount = currentParticle.mergedCount;
					upParticle.mergedProduct = 1;

					downParticle.parentData = currentParticle.ownedData;
					downParticle.hasNextEdge = false;
					downParticle.hasNextEdgeVector.push_back(false);
					downParticle.parentIndices.push_back(particleCounter);
					downParticle.weight = currentParticle.weight * inopProbability;
					downParticle.importanceDensity = currentParticle.importanceDensity * qTilde;
					downParticle.minCutSize = minCutSizeDown;
					downParticle.mergedCount = currentParticle.mergedCount;
					downParticle.mergedProduct = 1;

					newParticles.emplace_back(std::move(downParticle));
					newParticles.emplace_back(std::move(upParticle));
				}
			}
			//Mark as active edges whoose state is irrelevant at this point. 
			for(int particleCounter = 0; particleCounter < (int)newParticles.size(); particleCounter++)
			{
				approximateZeroVarianceWORMergeWithVarianceImpl::particle& currentParticle = newParticles[particleCounter];
				typedef boost::color_traits<boost::default_color_type> Color;
				std::fill(scratch.colorVector.begin(), scratch.colorVector.end(), Color::white());
				//compute connected components. 
				if(currentParticle.hasNextEdge)
				{
					currentParticle.parentData->capacity[2*edgeCounter] = currentParticle.parentData->capacity[2*edgeCounter+1] = HIGH_CAPACITY;
				}
				else currentParticle.parentData->capacity[2*edgeCounter] = currentParticle.parentData->capacity[2*edgeCounter+1] = 0;
				boost::connected_components_fixed(undirectedGraph, &(components[0]), &(scratch.colorVector[0]), componentsStack, &(currentParticle.parentData->capacity[0]));
				//Check which vertices are accessible from the source or sink, via edges with capacity HIGH_CAPACITY
				bool altered = false;
				std::vector<int>& downEdges = currentParticle.parentData->downEdges;
				std::vector<int> newDownEdges;
				for(int i = 0; i < (int)downEdges.size(); i++)
				{
					int edgeIndex = downEdges[i];
					int vertex1 = edges[edgeIndex].first;
					int vertex2 = edges[edgeIndex].second;
					if(components[vertex1] == components[vertex2])
					{
						if(!altered)
						{
							altered = true;
							currentParticle.ownedData.reset(new approximateZeroVarianceWORMergeWithVarianceImpl::particleData());
							currentParticle.ownedData->capacity = currentParticle.parentData->capacity;
							currentParticle.ownedData->residual = currentParticle.parentData->residual;
							if(currentParticle.hasNextEdge)
							{
								currentParticle.ownedData->capacity[2*edgeCounter] = currentParticle.ownedData->capacity[2*edgeCounter+1] = currentParticle.ownedData->residual[2*edgeCounter] = currentParticle.ownedData->residual[2*edgeCounter+1] = HIGH_CAPACITY;
							}
							else
							{
								currentParticle.ownedData->capacity[2*edgeCounter] = currentParticle.ownedData->capacity[2*edgeCounter+1] = currentParticle.ownedData->residual[2*edgeCounter] = currentParticle.ownedData->residual[2*edgeCounter+1] = 0;
							}
						}
						currentParticle.ownedData->capacity[2 * edgeIndex] = currentParticle.ownedData->capacity[2 * edgeIndex + 1] = HIGH_CAPACITY;
						currentParticle.ownedData->residual[2 * edgeIndex] = currentParticle.ownedData->residual[2 * edgeIndex + 1] = HIGH_CAPACITY;
					}
					else
					{
						newDownEdges.push_back(edgeIndex);
					}
				}
				if(!currentParticle.hasNextEdge)
				{
					int vertex1 = edges[edgeCounter].first;
					int vertex2 = edges[edgeCounter].second;
					if(components[vertex1] == components[vertex2])
					{
						currentParticle.hasNextEdge = true;
					}
				}
				if(altered) currentParticle.ownedData->downEdges = newDownEdges;
			}
			//Sort by state.
			std::sort(newParticles.begin(), newParticles.end(), [edgeCounter](const approximateZeroVarianceWORMergeWithVarianceImpl::particle& first, const approximateZeroVarianceWORMergeWithVarianceImpl::particle& second){ return first.order(second, edgeCounter);});
			std::size_t unitsAfterMerge = 0;
			//Merge particles
			{
				int particleCounter = 0;
				while(particleCounter < (int)newParticles.size())
				{
					mpfr_class additionalWeight = 0, additionalImportanceDensity = 0;
					int mergedCount = 1, additionalMergedProduct = 0;
					int particleCounter2 = particleCounter+1;
					approximateZeroVarianceWORMergeWithVarianceImpl::particle& currentParticle = newParticles[particleCounter];
					for(; particleCounter2 < (int)newParticles.size() && newParticles[particleCounter].matches(newParticles[particleCounter2], edgeCounter); particleCounter2++)
					{
						approximateZeroVarianceWORMergeWithVarianceImpl::particle& otherParticle = newParticles[particleCounter2];
						additionalWeight += otherParticle.weight;
						additionalImportanceDensity += otherParticle.importanceDensity;
						additionalMergedProduct += otherParticle.mergedProduct;
						mergedCount++;
						currentParticle.parentIndices.push_back(otherParticle.parentIndices[0]);
						currentParticle.hasNextEdgeVector.push_back(otherParticle.hasNextEdgeVector[0]);
						otherParticle.weight = 0;
						otherParticle.importanceDensity = 0;
						otherParticle.mergedCount = otherParticle.mergedProduct = 0;
					}
					currentParticle.weight += additionalWeight;
					currentParticle.importanceDensity += additionalImportanceDensity;
					currentParticle.mergedProduct += additionalMergedProduct;
					currentParticle.mergedCount = mergedCount;
					unitsAfterMerge++;
					particleCounter = particleCounter2;
				}
			}
			//std::size_t unitsAfterMerge = newParticles.size();
			samplingArgs.weights.clear();
			for(std::vector<approximateZeroVarianceWORMergeWithVarianceImpl::particle>::iterator i = newParticles.begin(); i != newParticles.end(); i++)
			{
				samplingArgs.weights.push_back(i->importanceDensity);
			}
			particles.clear();
			if(unitsAfterMerge <= n)
			{
				int counter = 0;
				for(int i = 0; i < (int)newParticles.size(); i++)
				{
					if(newParticles[i].weight != 0) 
					{
						particles.emplace_back(std::move(newParticles[i]));
						approximateZeroVarianceWORMergeWithVarianceImpl::particle& movedParticle = particles.back();
						if(!(movedParticle.ownedData))
						{
							movedParticle.ownedData.reset(new approximateZeroVarianceWORMergeWithVarianceImpl::particleData());
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
						approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraph::vertex_descriptor newVertex = boost::add_vertex(varianceEstimationGraph);
						approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraphVertex& vertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, newVertex);
						vertexInfo.samplingStage = edgeCounter;
						vertexInfo.indexWithinDesign = counter;
						vertexInfo.mergedCount = movedParticle.mergedCount;
						vertexInfo.mergedProduct = movedParticle.mergedProduct;
						for(int i = 0; i < (int)movedParticle.parentIndices.size(); i++)
						{
							boost::add_edge(graphVertices(edgeCounter-1, movedParticle.parentIndices[i]), newVertex, (bool)movedParticle.hasNextEdgeVector[i], varianceEstimationGraph);
						}
						graphVertices(edgeCounter, counter) = (int)newVertex;
						counter++;
					}
				}
				std::vector<::sampling::mpfr_class> inclusion(unitsAfterMerge, 1.0);
				allInclusionProbabilities.emplace_back(std::move(inclusion));
				boost::numeric::ublas::matrix<mpfr_class> secondOrder(unitsAfterMerge, unitsAfterMerge, 1.0);
				allSecondOrderInclusionProbabilities.emplace_back(std::move(secondOrder));
			}
			else
			{
				indices.clear();
				sampling::conditionalPoissonSequential(samplingArgs, args.randomSource);
				int counter = 0;
				for(std::vector<int>::iterator i = indices.begin(); i != indices.end(); i++)
				{
					particles.emplace_back(std::move(newParticles[*i]));
					approximateZeroVarianceWORMergeWithVarianceImpl::particle& movedParticle = particles.back();
					if(!(movedParticle.ownedData))
					{
						movedParticle.ownedData.reset(new approximateZeroVarianceWORMergeWithVarianceImpl::particleData());
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
					movedParticle.importanceDensity /= samplingArgs.inclusionProbabilities[*i];
					movedParticle.weight /= samplingArgs.inclusionProbabilities[*i];

					approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraph::vertex_descriptor newVertex = boost::add_vertex(varianceEstimationGraph);
					approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraphVertex& vertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, newVertex);
					vertexInfo.samplingStage = edgeCounter;
					vertexInfo.indexWithinDesign = *i;
					vertexInfo.mergedCount = movedParticle.mergedCount;
					vertexInfo.mergedProduct = movedParticle.mergedProduct;
					for(int i = 0; i < (int)movedParticle.parentIndices.size(); i++)
					{
						boost::add_edge(graphVertices(edgeCounter-1, movedParticle.parentIndices[i]), newVertex, (bool)movedParticle.hasNextEdgeVector[i], varianceEstimationGraph);
					}
					graphVertices(edgeCounter, counter) = (int)newVertex;
					counter++;
				}
				boost::numeric::ublas::matrix<mpfr_class> secondOrder;
				sampling::conditionalPoissonSecondOrderInclusionProbabilities(samplingArgs, samplingArgs.inclusionProbabilities, secondOrder);
				allSecondOrderInclusionProbabilities.emplace_back(std::move(secondOrder));
				allInclusionProbabilities.emplace_back(std::move(samplingArgs.inclusionProbabilities));
			}
		}
		//Put actual estimate values into the graph
		for(int i = 0; i < (int)particles.size(); i++)
		{
			if(graphVertices(nEdges - 1, i) > -1)
			{
				approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraphVertex& vertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, graphVertices(nEdges - 1, i));
				vertexInfo.accumulatedMean = 1;
			}
		}
		//Now for the covariance matrix estimations. 
		std::vector<boost::numeric::ublas::matrix<mpfr_class> > allCovariances(nEdges);
		for(int edgeCounter = nEdges - 2; edgeCounter >= 0; edgeCounter--)
		{
			int particleCounter = 0;
			const std::vector<mpfr_class>& currentInclusionProbabilities = allInclusionProbabilities[edgeCounter+1];
			while(particleCounter < (int)n && graphVertices(edgeCounter, particleCounter) > -1)
			{
				int currentVertex = graphVertices(edgeCounter, particleCounter);
				approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraphVertex& vertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, currentVertex);
				approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraph::out_edge_iterator current, end;
				boost::tie(current, end) = boost::out_edges(currentVertex, varianceEstimationGraph);
				for(; current != end; current++)
				{
					int targetVertex = boost::target(*current, varianceEstimationGraph);
					approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraphVertex& targetVertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, targetVertex);
					if(boost::get(boost::edge_name, varianceEstimationGraph, *current))
					{
						vertexInfo.accumulatedMean += targetVertexInfo.accumulatedMean * opProbability / currentInclusionProbabilities[targetVertexInfo.indexWithinDesign];
					}
					else
					{
						vertexInfo.accumulatedMean += targetVertexInfo.accumulatedMean * inopProbability / currentInclusionProbabilities[targetVertexInfo.indexWithinDesign];
					}
				}
				particleCounter++;
			}
		}
		for(std::vector<approximateZeroVarianceWORMergeWithVarianceImpl::particle>::iterator i = particles.begin(); i != particles.end(); i++)
		{
			args.estimate += i->weight;
		}
		//mpfr_class totalFromGraph = boost::get(boost::vertex_name, varianceEstimationGraph, 1).accumulatedMean*inopProbability + boost::get(boost::vertex_name, varianceEstimationGraph, 2).accumulatedMean*opProbability;
	}
}
