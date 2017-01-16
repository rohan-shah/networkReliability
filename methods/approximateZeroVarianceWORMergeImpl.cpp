#include "approximateZeroVarianceWORImpl.h"
#include <boost/random/random_number_generator.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include "sampford.h"
#include "depth_first_search_fixed.hpp"
namespace networkReliability
{
	namespace approximateZeroVarianceWORMergeImpl
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
		struct particle
		{
		public:
			boost::shared_ptr<particleData> parentData;
			boost::shared_ptr<particleData> ownedData;
			mpfr_class weight, importanceDensity;
			bool hasNextEdge;
			int minCutSize;
			particle()
			{}
			particle(particle&& other)
				: parentData(std::move(other.parentData)), ownedData(std::move(other.ownedData)), weight(std::move(other.weight)), importanceDensity(std::move(other.importanceDensity)), hasNextEdge(other.hasNextEdge), minCutSize(other.minCutSize)
			{}
			particle& operator=(particle&& other)
			{
				parentData.swap(other.parentData);
				ownedData.swap(other.ownedData);
				weight = other.weight;
				importanceDensity = other.importanceDensity;
				hasNextEdge = other.hasNextEdge;
				minCutSize = other.minCutSize;
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
	int getMinCut(std::vector<int>& capacity, std::vector<int>& residual, const context::internalDirectedGraph& graph, const context::internalGraph& undirectedGraph, const std::vector<int>& interestVertices, approximateZeroVarianceWORMergeImpl::approximateZeroVarianceScratch& scratch)
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
	void approximateZeroVarianceWORMerge(approximateZeroVarianceWORArgs& args)
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

		boost::detail::depth_first_visit_fixed_impl_helper<context::internalGraph>::stackType fixedSearchStack;
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
		approximateZeroVarianceWORMergeImpl::approximateZeroVarianceScratch scratch;
		scratch.colorVector.resize(nVertices);

		//Initialise with the two initial choices - The first edge can be up or down. 
		std::vector<mpfr_class>& newImportanceDensity = samplingArgs.weights;
		std::vector<mpfr_class>& rescaledWeights = samplingArgs.rescaledWeights;

		std::vector<approximateZeroVarianceWORMergeImpl::particle> particles, newParticles;
		{
			approximateZeroVarianceWORMergeImpl::particle missingEdgeParticle;
			missingEdgeParticle.ownedData.reset(new approximateZeroVarianceWORMergeImpl::particleData());
			missingEdgeParticle.hasNextEdge = false;
			missingEdgeParticle.ownedData->capacity.resize(2*nEdges, 1);
			missingEdgeParticle.ownedData->residual.resize(2*nEdges, 1);
			missingEdgeParticle.ownedData->capacity[0] = missingEdgeParticle.ownedData->capacity[1] = 0;
			missingEdgeParticle.ownedData->residual[0] = missingEdgeParticle.ownedData->residual[1] = 0;
			missingEdgeParticle.minCutSize = getMinCut(missingEdgeParticle.ownedData->capacity, missingEdgeParticle.ownedData->residual, graph, undirectedGraph, interestVertices, scratch);
			missingEdgeParticle.weight = inopProbability;
			missingEdgeParticle.ownedData->downEdges.push_back(0);

			approximateZeroVarianceWORMergeImpl::particle withEdgeParticle;
			withEdgeParticle.ownedData.reset(new approximateZeroVarianceWORMergeImpl::particleData());
			withEdgeParticle.hasNextEdge = true;
			withEdgeParticle.ownedData->capacity.resize(2*nEdges, 1);
			withEdgeParticle.ownedData->residual.resize(2*nEdges, 1);
			withEdgeParticle.ownedData->capacity[0] = withEdgeParticle.ownedData->capacity[1] = HIGH_CAPACITY;
			withEdgeParticle.ownedData->residual[0] = withEdgeParticle.ownedData->residual[1] = HIGH_CAPACITY;
			withEdgeParticle.minCutSize = getMinCut(withEdgeParticle.ownedData->capacity, withEdgeParticle.ownedData->residual, graph, undirectedGraph, interestVertices, scratch);
			withEdgeParticle.weight = opProbability;

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
			}
			else
			{
				missingEdgeParticle.importanceDensity = 1;
				particles.emplace_back(std::move(missingEdgeParticle));
			}
		}
		args.estimate = 0;
		for(int edgeCounter = 1; edgeCounter < (int)nEdges; edgeCounter++)
		{
			newParticles.clear();
			for(int particleCounter = 0; particleCounter < (int)particles.size(); particleCounter++)
			{
				approximateZeroVarianceWORMergeImpl::particle& currentParticle = particles[particleCounter];
				if(currentParticle.minCutSize == 0)
				{
					args.estimate += currentParticle.weight;
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
					approximateZeroVarianceWORMergeImpl::particle newParticle;
					newParticle.parentData = particles[particleCounter].ownedData;
					newParticle.hasNextEdge = false;
					newParticle.weight = currentParticle.weight * inopProbability;
					newParticle.importanceDensity = currentParticle.importanceDensity;
					newParticle.minCutSize = minCutSizeDown;
					newParticles.emplace_back(std::move(newParticle));
				}
				else
				{
					mpfr_class minCutDownProb = cachedInopPowers[minCutSizeDown];
					mpfr_class minCutUpProb = cachedInopPowers[minCutSizeUp];
					mpfr_class qTilde = inopProbability * minCutDownProb;
					qTilde = qTilde / (qTilde  + opProbability * minCutUpProb);
					
					approximateZeroVarianceWORMergeImpl::particle upParticle, downParticle;

					upParticle.parentData = currentParticle.ownedData;
					upParticle.hasNextEdge = true;
					upParticle.weight = currentParticle.weight * opProbability;
					upParticle.importanceDensity = currentParticle.importanceDensity * (1 - qTilde);
					upParticle.minCutSize = minCutSizeUp;

					downParticle.parentData = currentParticle.ownedData;
					downParticle.hasNextEdge = false;
					downParticle.weight = currentParticle.weight * inopProbability;
					downParticle.importanceDensity = currentParticle.importanceDensity * qTilde;
					downParticle.minCutSize = minCutSizeDown;

					newParticles.emplace_back(std::move(downParticle));
					newParticles.emplace_back(std::move(upParticle));
				}
			}
			//Mark as active edges whoose state is irrelevant at this point. 
			for(int particleCounter = 0; particleCounter < (int)newParticles.size(); particleCounter++)
			{
				approximateZeroVarianceWORMergeImpl::particle& currentParticle = newParticles[particleCounter];
				typedef boost::color_traits<boost::default_color_type> Color;
				std::fill(scratch.colorVector.begin(), scratch.colorVector.end(), Color::white());
				boost::default_dfs_visitor visitor;
				//Check which vertices are accessible from the source or sink, via edges with capacity HIGH_CAPACITY
				if(currentParticle.ownedData)
				{
					throw std::runtime_error("This should be impossible");
				}
				else
				{
					int copied = currentParticle.parentData->capacity[2*edgeCounter];
					if(currentParticle.hasNextEdge)
					{
						currentParticle.parentData->capacity[2*edgeCounter] = currentParticle.parentData->capacity[2*edgeCounter+1] = HIGH_CAPACITY;
					}
					else currentParticle.parentData->capacity[2*edgeCounter] = currentParticle.parentData->capacity[2*edgeCounter+1] = 0;
					boost::detail::depth_first_visit_fixed_impl(undirectedGraph, interestVertices[0], visitor, &(scratch.colorVector[0]), fixedSearchStack, &(currentParticle.parentData->capacity[0]), boost::detail::nontruth2());
					currentParticle.parentData->capacity[2*edgeCounter] = currentParticle.parentData->capacity[2*edgeCounter+1] = copied;
				}
				bool altered = false;
				std::vector<int>& downEdges = currentParticle.parentData->downEdges;
				for(int i = 0; i < (int)downEdges.size(); i++)
				{
					int edgeIndex = downEdges[i];
					int vertex1 = edges[edgeIndex].first;
					int vertex2 = edges[edgeIndex].second;
					if(scratch.colorVector[vertex1] == Color::black() && scratch.colorVector[vertex2] == Color::black())
					{
						if(!altered)
						{
							altered = true;
							currentParticle.ownedData.reset(new approximateZeroVarianceWORMergeImpl::particleData());
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
							std::copy(currentParticle.parentData->downEdges.begin(), currentParticle.parentData->downEdges.begin() + i, currentParticle.ownedData->downEdges.begin());
						}
						currentParticle.ownedData->capacity[2 * edgeIndex] = currentParticle.ownedData->capacity[2 * edgeIndex + 1] = HIGH_CAPACITY;
						currentParticle.ownedData->residual[2 * edgeIndex] = currentParticle.ownedData->residual[2 * edgeIndex + 1] = HIGH_CAPACITY;
					}
					else if(altered)
					{
						currentParticle.ownedData->downEdges.push_back(edgeIndex);
					}
				}


				std::fill(scratch.colorVector.begin(), scratch.colorVector.end(), Color::white());
				if(currentParticle.ownedData)
				{
					boost::detail::depth_first_visit_fixed_impl(undirectedGraph, interestVertices[1], visitor, &(scratch.colorVector[0]), fixedSearchStack, &(currentParticle.ownedData->capacity[0]), boost::detail::nontruth2());
				}
				else
				{
					int copied = currentParticle.parentData->capacity[2*edgeCounter];
					if(currentParticle.hasNextEdge)
					{
						currentParticle.parentData->capacity[2*edgeCounter] = currentParticle.parentData->capacity[2*edgeCounter+1] = HIGH_CAPACITY;
					}
					else currentParticle.parentData->capacity[2*edgeCounter] = currentParticle.parentData->capacity[2*edgeCounter+1] = 0;
					boost::detail::depth_first_visit_fixed_impl(undirectedGraph, interestVertices[1], visitor, &(scratch.colorVector[0]), fixedSearchStack, &(currentParticle.parentData->capacity[0]), boost::detail::nontruth2());
					currentParticle.parentData->capacity[2*edgeCounter] = currentParticle.parentData->capacity[2*edgeCounter+1] = copied;
				}
				for(int i = 0; i < (int)downEdges.size(); i++)
				{
					int edgeIndex = downEdges[i];
					int vertex1 = edges[edgeIndex].first;
					int vertex2 = edges[edgeIndex].second;
					if(scratch.colorVector[vertex1] == Color::black() && scratch.colorVector[vertex2] == Color::black())
					{
						if(!altered)
						{
							altered = true;
							currentParticle.ownedData.reset(new approximateZeroVarianceWORMergeImpl::particleData());
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
							std::copy(currentParticle.parentData->downEdges.begin(), currentParticle.parentData->downEdges.begin() + i, currentParticle.ownedData->downEdges.begin());
						}
						currentParticle.ownedData->capacity[2 * edgeIndex] = currentParticle.ownedData->capacity[2 * edgeIndex + 1] = HIGH_CAPACITY;
						currentParticle.ownedData->residual[2 * edgeIndex] = currentParticle.ownedData->residual[2 * edgeIndex + 1] = HIGH_CAPACITY;
					}
					else if(altered)
					{
						currentParticle.ownedData->downEdges.push_back(edgeIndex);
					}
				}
			}
			//Sort by state.
			std::sort(newParticles.begin(), newParticles.end(), [edgeCounter](const approximateZeroVarianceWORMergeImpl::particle& first, const approximateZeroVarianceWORMergeImpl::particle& second){ return first.order(second, edgeCounter);});
			std::size_t unitsAfterMerge = 0;
			//Merge particles
			{
				int particleCounter = 0;
				while(particleCounter < (int)newParticles.size())
				{
					mpfr_class additionalWeight = 0, additionalImportanceDensity = 0;
					int particleCounter2 = particleCounter+1;
					for(; particleCounter2 < (int)newParticles.size() && newParticles[particleCounter].matches(newParticles[particleCounter2], edgeCounter); particleCounter2++)
					{
						additionalWeight += newParticles[particleCounter2].weight;
						additionalImportanceDensity += newParticles[particleCounter2].importanceDensity;
						newParticles[particleCounter2].weight = 0;
						newParticles[particleCounter2].importanceDensity = 0;
					}
					newParticles[particleCounter].weight += additionalWeight;
					newParticles[particleCounter].importanceDensity += additionalImportanceDensity;
					unitsAfterMerge++;
					particleCounter = particleCounter2;
				}
			}
			//std::size_t unitsAfterMerge = newParticles.size();
			newImportanceDensity.clear();
			for(std::vector<approximateZeroVarianceWORMergeImpl::particle>::iterator i = newParticles.begin(); i != newParticles.end(); i++)
			{
				newImportanceDensity.push_back(i->importanceDensity);
			}
			particles.clear();
			if(unitsAfterMerge <= n)
			{
				for(int i = 0; i < (int)newParticles.size(); i++)
				{
					if(newParticles[i].weight != 0) 
					{
						particles.emplace_back(std::move(newParticles[i]));
						approximateZeroVarianceWORMergeImpl::particle& movedParticle = particles.back();
						if(!(movedParticle.ownedData))
						{
							movedParticle.ownedData.reset(new approximateZeroVarianceWORMergeImpl::particleData());
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
							}
						}
					}
				}
			}
			else
			{
				indices.clear();
				sampling::sampfordFromParetoNaive(samplingArgs, args.randomSource);
				for(std::vector<int>::iterator i = indices.begin(); i != indices.end(); i++)
				{
					particles.emplace_back(std::move(newParticles[*i]));
					approximateZeroVarianceWORMergeImpl::particle& movedParticle = particles.back();
					if(!(movedParticle.ownedData))
					{
						movedParticle.ownedData.reset(new approximateZeroVarianceWORMergeImpl::particleData());
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
						}
					}
					movedParticle.importanceDensity /= rescaledWeights[*i];
					movedParticle.weight /= rescaledWeights[*i];
				}
			}
		}
		for(std::vector<approximateZeroVarianceWORMergeImpl::particle>::iterator i = particles.begin(); i != particles.end(); i++)
		{
			args.estimate += i->weight;
		}
	}
}
