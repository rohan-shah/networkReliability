#include "approximateZeroVarianceWORMergeWithVarianceImpl.h"
#include <boost/random/random_number_generator.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include "conditionalPoissonSequential.h"
#include "connected_components_fixed.hpp"
#include <boost/graph/graphml.hpp>
#include <boost/graph/graphviz.hpp>
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
				: indexWithinDesign(-1), samplingStage(-1), mergedCount(-1), mergedProduct(-1), accumulatedMean(0), trueDensity(0)//, V(0)
			{}
			int indexWithinDesign;
			int samplingStage;
			int mergedCount;
			int mergedProduct;
			int indexWithinSelected;
			mutable ::sampling::mpfr_class accumulatedMean, trueDensity;//, V;
		};
		typedef boost::adjacency_list<boost::listS, boost::vecS, boost::directedS, boost::property<boost::vertex_name_t, varianceGraphVertex>, boost::property<boost::edge_name_t, bool> > varianceGraph;
		struct particle
		{
		public:
			boost::shared_ptr<particleData> parentData;
			boost::shared_ptr<particleData> ownedData;
			mpfr_class weight, importanceDensity, trueDensity;
			bool hasNextEdge;
			int minCutSize;
			int mergedCount;
			int mergedProduct;
			std::vector<int> parentIndices;
			std::vector<bool> hasNextEdgeVector;
			particle()
			{}
			particle(particle&& other)
				: parentData(std::move(other.parentData)), ownedData(std::move(other.ownedData)), weight(std::move(other.weight)), importanceDensity(std::move(other.importanceDensity)), trueDensity(std::move(other.trueDensity)), hasNextEdge(other.hasNextEdge), minCutSize(other.minCutSize), mergedCount(other.mergedCount), mergedProduct(other.mergedProduct), parentIndices(std::move(other.parentIndices)), hasNextEdgeVector(std::move(other.hasNextEdgeVector))
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
				trueDensity = other.trueDensity;
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
/*	class vertexPropertyWriter
	{
	public:
		vertexPropertyWriter(const approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraph& g)
			:g(g)
		{}
		template<typename vertexDesc> void operator()(std::ostream& out, const vertexDesc& v)
		{
			const approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraphVertex& vertexInfo = boost::get(boost::vertex_name, g, v);
			out << std::endl << "[trueDensity=\"" << vertexInfo.trueDensity.convert_to<double>() << "\"";
			out << ",accumulatedMean=\"" << vertexInfo.accumulatedMean.convert_to<double>() << "\"";
			out << ",samplingStage=\"" << vertexInfo.samplingStage << "\"";
			out << ",indexWithinDesign=\"" << vertexInfo.indexWithinDesign << "\"";
			out << ",V=\"" << vertexInfo.V.convert_to<double>() << "\"";
			out << ",indexWithinSelected=\"" << vertexInfo.indexWithinSelected <<"\"]" << std::endl;
		}
		const approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraph& g;
	};
	class edgePropertyWriter
	{
	public:
		edgePropertyWriter(const approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraph& g)
			:g(g)
		{}
		template<typename edgeDesc> void operator()(std::ostream& out, const edgeDesc& e)
		{
			bool edgeInfo = boost::get(boost::edge_name, g, e);
			out << " [isPresent=\"" << edgeInfo <<"\"]" << std::endl;
		}
		const approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraph& g;
	};*/
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
			missingEdgeParticle.trueDensity = missingEdgeParticle.weight = inopProbability;
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
			withEdgeParticle.trueDensity = withEdgeParticle.weight = opProbability;
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
				firstVertexInfo.indexWithinDesign = firstVertexInfo.indexWithinSelected = 0;
				firstVertexInfo.mergedCount = firstVertexInfo.mergedProduct = 1;
				firstVertexInfo.trueDensity = inopProbability;
				boost::add_edge(rootVertex, firstVertex, varianceEstimationGraph);
				graphVertices(0, 0) = (int)firstVertex;

				approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraph::vertex_descriptor secondVertex = boost::add_vertex(varianceEstimationGraph);
				approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraphVertex& secondVertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, secondVertex);
				secondVertexInfo.samplingStage = 0;
				secondVertexInfo.indexWithinDesign = secondVertexInfo.indexWithinSelected = 1;
				secondVertexInfo.mergedCount = secondVertexInfo.mergedProduct = 1;
				secondVertexInfo.trueDensity = opProbability;
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
				vertexInfo.indexWithinDesign = vertexInfo.indexWithinSelected = 0;
				vertexInfo.mergedCount = vertexInfo.mergedProduct = 1;
				vertexInfo.trueDensity = inopProbability;
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
					vertexInfo.mergedCount = currentParticle.mergedCount;
					vertexInfo.mergedProduct = currentParticle.mergedProduct;
					vertexInfo.trueDensity = currentParticle.trueDensity;
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
					newParticle.trueDensity = currentParticle.trueDensity * inopProbability;
					newParticle.importanceDensity = currentParticle.importanceDensity;
					newParticle.minCutSize = minCutSizeDown;
					newParticle.parentIndices.push_back(particleCounter);
					newParticles.emplace_back(std::move(newParticle));
					newParticle.mergedCount = 1;
					newParticle.mergedProduct = currentParticle.mergedProduct;
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
					upParticle.trueDensity = currentParticle.trueDensity * opProbability;
					upParticle.importanceDensity = currentParticle.importanceDensity * (1 - qTilde);
					upParticle.minCutSize = minCutSizeUp;
					upParticle.mergedCount = 1;
					upParticle.mergedProduct = currentParticle.mergedCount;

					downParticle.parentData = currentParticle.ownedData;
					downParticle.hasNextEdge = false;
					downParticle.hasNextEdgeVector.push_back(false);
					downParticle.parentIndices.push_back(particleCounter);
					downParticle.weight = currentParticle.weight * inopProbability;
					downParticle.trueDensity = currentParticle.trueDensity * inopProbability;
					downParticle.importanceDensity = currentParticle.importanceDensity * qTilde;
					downParticle.minCutSize = minCutSizeDown;
					downParticle.mergedCount = 1;
					downParticle.mergedProduct = currentParticle.mergedProduct;

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
					mpfr_class additionalWeight = 0, additionalImportanceDensity = 0, additionalTrueDensity = 0;
					int mergedCount = 1, additionalMergedProduct = 0;
					int particleCounter2 = particleCounter+1;
					approximateZeroVarianceWORMergeWithVarianceImpl::particle& currentParticle = newParticles[particleCounter];
					for(; particleCounter2 < (int)newParticles.size() && newParticles[particleCounter].matches(newParticles[particleCounter2], edgeCounter); particleCounter2++)
					{
						approximateZeroVarianceWORMergeWithVarianceImpl::particle& otherParticle = newParticles[particleCounter2];
						additionalWeight += otherParticle.weight;
						additionalImportanceDensity += otherParticle.importanceDensity;
						additionalTrueDensity += otherParticle.trueDensity;
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
					currentParticle.trueDensity += additionalTrueDensity;
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
						vertexInfo.indexWithinDesign = vertexInfo.indexWithinSelected = counter;
						vertexInfo.mergedCount = movedParticle.mergedCount;
						vertexInfo.mergedProduct = movedParticle.mergedProduct;
						vertexInfo.trueDensity = movedParticle.trueDensity;
						for(int i = 0; i < (int)movedParticle.parentIndices.size(); i++)
						{
							boost::add_edge(graphVertices(edgeCounter-1, movedParticle.parentIndices[i]), newVertex, (bool)movedParticle.hasNextEdgeVector[i], varianceEstimationGraph);
						}
						movedParticle.parentIndices.clear();
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
				bool hasTooLarge = true;
				while(hasTooLarge)
				{
					hasTooLarge = false;
					sampling::conditionalPoissonSequential(samplingArgs, args.randomSource);
					int nDeterministic = 0;
					for(int i = 0; i < (int) samplingArgs.inclusionProbabilities.size(); i++)
					{
						if(samplingArgs.inclusionProbabilities[i] > 0.9999 && samplingArgs.inclusionProbabilities[i] != 1)
						{
							hasTooLarge = true;
							samplingArgs.weights[i] = 1.0;
						}
						if(samplingArgs.weights[i] == 1) nDeterministic++;
					}
					samplingArgs.n = std::max(nDeterministic+1, (int)n);
				}
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
					vertexInfo.indexWithinSelected = (int)std::distance(indices.begin(), i);
					vertexInfo.mergedCount = movedParticle.mergedCount;
					vertexInfo.mergedProduct = movedParticle.mergedProduct;
					vertexInfo.trueDensity = movedParticle.trueDensity;
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
		//Put in values of 1 for the accumulated means at the end
		for(int i = 0; i < (int)particles.size(); i++)
		{
			if(graphVertices(nEdges - 1, i) > -1)
			{
				approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraphVertex& vertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, graphVertices(nEdges - 1, i));
				vertexInfo.accumulatedMean = 1;
			}
		}
		//Now for the covariance matrix estimations. 
		std::vector<boost::numeric::ublas::matrix<mpfr_class> > allCovariances(nEdges, boost::numeric::ublas::matrix<mpfr_class>(n, n, 0));
		for(int edgeCounter = nEdges - 2; edgeCounter >= 0; edgeCounter--)
		{
			int particleCounter = 0;
			const std::vector<mpfr_class>& currentInclusionProbabilities = allInclusionProbabilities[edgeCounter+1];
			const boost::numeric::ublas::matrix<mpfr_class>& currentSecondOrderInclusionProbabilities = allSecondOrderInclusionProbabilities[edgeCounter+1];
			boost::numeric::ublas::matrix<mpfr_class>& currentCovariance = allCovariances[edgeCounter];
			const boost::numeric::ublas::matrix<mpfr_class>& previousCovariance = allCovariances[edgeCounter+1];
			//First estimate the accumulated means
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
						vertexInfo.accumulatedMean += (opProbability * targetVertexInfo.accumulatedMean / currentInclusionProbabilities[targetVertexInfo.indexWithinDesign]);
					}
					else
					{
						vertexInfo.accumulatedMean += (inopProbability * targetVertexInfo.accumulatedMean / currentInclusionProbabilities[targetVertexInfo.indexWithinDesign]);
					}
				}
				particleCounter++;
			}
			//now actually do the covariance estimation. 
			int particleCounter1 = 0;
			while(particleCounter1 < (int)n && graphVertices(edgeCounter, particleCounter1) > -1)
			{
				int particleCounter2 = 0;
				while(particleCounter2 < (int)n && graphVertices(edgeCounter, particleCounter2) > -1)
				{
					int graphVertex1 = graphVertices(edgeCounter, particleCounter1);
					approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraphVertex& graphVertex1Info = boost::get(boost::vertex_name, varianceEstimationGraph, graphVertex1);
					int graphVertex2 = graphVertices(edgeCounter, particleCounter2);
					approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraphVertex& graphVertex2Info = boost::get(boost::vertex_name, varianceEstimationGraph, graphVertex2);

					mpfr_class& currentCovarianceValue = currentCovariance(particleCounter1, particleCounter2);
				
					approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraph::out_edge_iterator current1, end1, current2, end2;
					boost::tie(current1, end1) = boost::out_edges(graphVertex1, varianceEstimationGraph);
					for(; current1 != end1; current1++)
					{
						int targetVertex1 = boost::target(*current1, varianceEstimationGraph);
						approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraphVertex& targetVertexInfo1 = boost::get(boost::vertex_name, varianceEstimationGraph, targetVertex1);
						mpfr_class multiple1;
						if(boost::get(boost::edge_name, varianceEstimationGraph, *current1)) multiple1 = opProbability;
						else multiple1 = inopProbability;

						boost::tie(current2, end2) = boost::out_edges(graphVertex2, varianceEstimationGraph);
						for(; current2 != end2; current2++)
						{
							int targetVertex2 = boost::target(*current2, varianceEstimationGraph);
							approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraphVertex& targetVertexInfo2 = boost::get(boost::vertex_name, varianceEstimationGraph, targetVertex2);
							mpfr_class inclusionProduct = currentInclusionProbabilities[targetVertexInfo1.indexWithinDesign] * currentInclusionProbabilities[targetVertexInfo2.indexWithinDesign];
							mpfr_class multiple2;
							if(boost::get(boost::edge_name, varianceEstimationGraph, *current2)) multiple2 = opProbability;
							else multiple2 = inopProbability;
							if(targetVertex1 == targetVertex2)
							{
								currentCovarianceValue += ((currentInclusionProbabilities[targetVertexInfo1.indexWithinDesign] - inclusionProduct) * targetVertexInfo1.accumulatedMean * targetVertexInfo2.accumulatedMean * graphVertex1Info.trueDensity * graphVertex2Info.trueDensity * multiple1 * multiple2 / (currentInclusionProbabilities[targetVertexInfo1.indexWithinDesign]*inclusionProduct)) + (previousCovariance(targetVertexInfo1.indexWithinSelected, targetVertexInfo2.indexWithinSelected) / (currentInclusionProbabilities[targetVertexInfo1.indexWithinDesign])) * (multiple1 * graphVertex1Info.trueDensity / targetVertexInfo1.trueDensity) * (multiple2 * graphVertex2Info.trueDensity / targetVertexInfo2.trueDensity);
							}
							else
							{
								currentCovarianceValue += ((currentSecondOrderInclusionProbabilities(targetVertexInfo1.indexWithinDesign, targetVertexInfo2.indexWithinDesign) - inclusionProduct) * targetVertexInfo1.accumulatedMean * targetVertexInfo2.accumulatedMean * graphVertex1Info.trueDensity * graphVertex2Info.trueDensity * multiple1 * multiple2 / (currentSecondOrderInclusionProbabilities(targetVertexInfo1.indexWithinDesign, targetVertexInfo2.indexWithinDesign) * inclusionProduct)) + (previousCovariance(targetVertexInfo1.indexWithinSelected, targetVertexInfo2.indexWithinSelected) / (currentSecondOrderInclusionProbabilities(targetVertexInfo1.indexWithinDesign, targetVertexInfo2.indexWithinDesign))) * (multiple1 * graphVertex1Info.trueDensity / targetVertexInfo1.trueDensity) * (multiple2 * graphVertex2Info.trueDensity / targetVertexInfo2.trueDensity);
							}
						}
					}
					/*if(particleCounter1 == particleCounter2)
					{
						graphVertex1Info.V = currentCovarianceValue;
					}*/
					particleCounter2++;
				}
				particleCounter1++;
			}
			//If a variance is zero, the covariances *have* to be zero. 
			particleCounter1 = 0;
			while(particleCounter1 < (int)n && graphVertices(edgeCounter, particleCounter1) > -1)
			{
				if(currentCovariance(particleCounter1, particleCounter1) == 0)
				{
					int particleCounter2 = 0;
					while(particleCounter2 < (int)n && graphVertices(edgeCounter, particleCounter2) > -1)
					{
						currentCovariance(particleCounter1, particleCounter2) = currentCovariance(particleCounter2, particleCounter1) = 0;
						particleCounter2++;
					}
				}
				particleCounter1++;
			}
		}
		for(std::vector<approximateZeroVarianceWORMergeWithVarianceImpl::particle>::iterator i = particles.begin(); i != particles.end(); i++)
		{
			args.estimate += i->weight;
		}
		mpfr_class totalFromGraph = boost::get(boost::vertex_name, varianceEstimationGraph, 1).accumulatedMean*inopProbability + boost::get(boost::vertex_name, varianceEstimationGraph, 2).accumulatedMean * opProbability;
		mpfr_class totalCovarianceFromGraph = allCovariances[0](0, 0) + allCovariances[0](1, 0) + allCovariances[0](0, 1) + allCovariances[0](1, 1);
		args.varianceEstimate = totalCovarianceFromGraph;
		mpfr_class relative1 = (totalFromGraph - args.estimate)/totalFromGraph;
		mpfr_class relative2 = (totalFromGraph - args.estimate)/args.estimate;
		if(std::fabs(relative1.convert_to<double>()) > 1e-3 || std::fabs(relative2.convert_to<double>()) > 1e-3)
		{
			throw std::runtime_error("Internal error");
		}
/*		{
			std::ofstream file("graph.graphml");
			boost::dynamic_properties dp;
			dp.property("samplingStage", boost::make_transform_value_property_map(boost::bind<int&>(&approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraphVertex::samplingStage, _1), boost::get(boost::vertex_name_t(), varianceEstimationGraph)));
			dp.property("indexWithinDesign", boost::make_transform_value_property_map(boost::bind<int&>(&approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraphVertex::indexWithinDesign, _1), boost::get(boost::vertex_name_t(), varianceEstimationGraph)));
			dp.property("indexWithinSelected", boost::make_transform_value_property_map(boost::bind<int&>(&approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraphVertex::indexWithinSelected, _1), boost::get(boost::vertex_name_t(), varianceEstimationGraph)));
			auto vertexToAccumulatedMean = [](approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraphVertex& x)
				{
					return x.accumulatedMean.convert_to<double>();
				};
			dp.property("accumulatedMean", boost::make_transform_value_property_map(vertexToAccumulatedMean, boost::get(boost::vertex_name_t(), varianceEstimationGraph)));
			auto vertexToTrueDensity = [](approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraphVertex& x)
				{
					return x.trueDensity.convert_to<double>();
				};
			dp.property("trueDensity", boost::make_transform_value_property_map(vertexToTrueDensity, boost::get(boost::vertex_name_t(), varianceEstimationGraph)));
			auto vertexToV = [](approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraphVertex& x)
				{
					return x.V.convert_to<double>();
				};
			dp.property("V", boost::make_transform_value_property_map(vertexToV, boost::get(boost::vertex_name_t(), varianceEstimationGraph)));
			boost::write_graphml(file, varianceEstimationGraph, dp);
		}
		{
			std::ofstream file("graph.dot");
			vertexPropertyWriter vp(varianceEstimationGraph);
			edgePropertyWriter ep(varianceEstimationGraph);
			boost::write_graphviz(file, varianceEstimationGraph, vp, ep);
		}*/
	}
}
