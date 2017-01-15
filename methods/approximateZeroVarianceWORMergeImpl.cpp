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
		struct particle
		{
		public:
			std::vector<int> capacity, residual, downEdges;
			mpfr_class weight, importanceDensity;
			int minCutSize;
			particle()
			{}
			particle(particle&& other)
				: capacity(std::move(other.capacity)), residual(std::move(other.residual)), downEdges(std::move(other.downEdges)), weight(std::move(other.weight)), importanceDensity(std::move(other.importanceDensity)), minCutSize(other.minCutSize)
			{}
			particle& operator=(particle&& other)
			{
				capacity.swap(other.capacity);
				residual.swap(other.residual);
				downEdges.swap(other.downEdges);
				weight = other.weight;
				importanceDensity = other.importanceDensity;
				minCutSize = other.minCutSize;
				return *this;
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
			missingEdgeParticle.capacity.resize(2*nEdges, 1);
			missingEdgeParticle.residual.resize(2*nEdges, 1);
			missingEdgeParticle.capacity[0] = missingEdgeParticle.capacity[1] = 0;
			missingEdgeParticle.residual[0] = missingEdgeParticle.residual[1] = 0;
			missingEdgeParticle.minCutSize = getMinCut(missingEdgeParticle.capacity, missingEdgeParticle.residual, graph, undirectedGraph, interestVertices, scratch);
			missingEdgeParticle.weight = inopProbability;
			missingEdgeParticle.downEdges.push_back(0);

			approximateZeroVarianceWORMergeImpl::particle withEdgeParticle;
			withEdgeParticle.capacity.resize(2*nEdges, 1);
			withEdgeParticle.residual.resize(2*nEdges, 1);
			withEdgeParticle.capacity[0] = withEdgeParticle.capacity[1] = HIGH_CAPACITY;
			withEdgeParticle.residual[0] = withEdgeParticle.residual[1] = HIGH_CAPACITY;
			withEdgeParticle.minCutSize = getMinCut(withEdgeParticle.capacity, withEdgeParticle.residual, graph, undirectedGraph, interestVertices, scratch);
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
				if(particles[particleCounter].minCutSize == 0)
				{
					args.estimate += particles[particleCounter].weight;
					continue;
				}
				particles[particleCounter].capacity[2*edgeCounter] = particles[particleCounter].capacity[2*edgeCounter + 1] = 0;
				particles[particleCounter].residual[2*edgeCounter] = particles[particleCounter].residual[2*edgeCounter + 1] = 0;
				int minCutSizeDown = getMinCut(particles[particleCounter].capacity, particles[particleCounter].residual, graph, undirectedGraph, interestVertices, scratch);

				particles[particleCounter].capacity[2*edgeCounter] = particles[particleCounter].capacity[2*edgeCounter + 1] = HIGH_CAPACITY;
				particles[particleCounter].residual[2*edgeCounter] = particles[particleCounter].residual[2*edgeCounter + 1] = HIGH_CAPACITY;
				int minCutSizeUp = getMinCut(particles[particleCounter].capacity, particles[particleCounter].residual, graph, undirectedGraph, interestVertices, scratch);
				if(minCutSizeUp >= HIGH_CAPACITY)
				{
					approximateZeroVarianceWORMergeImpl::particle newParticle;
					newParticle.capacity = particles[particleCounter].capacity;
					newParticle.capacity[2*edgeCounter] = newParticle.capacity[2*edgeCounter + 1] = 0;

					newParticle.residual = particles[particleCounter].residual;
					newParticle.residual[2*edgeCounter] = newParticle.residual[2*edgeCounter + 1] = 0;

					newParticle.weight = particles[particleCounter].weight * inopProbability;
					newParticle.importanceDensity = particles[particleCounter].importanceDensity;
					newParticle.minCutSize = minCutSizeDown;

					newParticle.downEdges = particles[particleCounter].downEdges;
					newParticle.downEdges.push_back(edgeCounter);
					newParticles.emplace_back(std::move(newParticle));
				}
				else
				{
					mpfr_class minCutDownProb = cachedInopPowers[minCutSizeDown];
					mpfr_class minCutUpProb = cachedInopPowers[minCutSizeUp];
					mpfr_class qTilde = inopProbability * minCutDownProb;
					qTilde = qTilde / (qTilde  + opProbability * minCutUpProb);
					
					approximateZeroVarianceWORMergeImpl::particle upParticle, downParticle;

					upParticle.capacity = particles[particleCounter].capacity;
					upParticle.capacity[2*edgeCounter] = upParticle.capacity[2*edgeCounter + 1] = HIGH_CAPACITY;
					upParticle.residual = particles[particleCounter].residual;
					upParticle.residual[2*edgeCounter] = upParticle.residual[2*edgeCounter + 1] = HIGH_CAPACITY;
					upParticle.weight = particles[particleCounter].weight * opProbability;
					upParticle.importanceDensity = particles[particleCounter].importanceDensity * (1 - qTilde);
					upParticle.minCutSize = minCutSizeUp;
					upParticle.downEdges = particles[particleCounter].downEdges;


					downParticle.capacity = particles[particleCounter].capacity;
					downParticle.capacity[2*edgeCounter] = downParticle.capacity[2*edgeCounter + 1] = 0;
					downParticle.residual = particles[particleCounter].residual;
					downParticle.residual[2*edgeCounter] = downParticle.residual[2*edgeCounter + 1] = 0;
					downParticle.weight = particles[particleCounter].weight * inopProbability;
					downParticle.importanceDensity = particles[particleCounter].importanceDensity * qTilde;
					downParticle.minCutSize = minCutSizeDown;
					downParticle.downEdges = particles[particleCounter].downEdges;
					downParticle.downEdges.push_back(edgeCounter);

					newParticles.emplace_back(std::move(downParticle));
					newParticles.emplace_back(std::move(upParticle));
				}
			}
			//Mark as active edges whoose state is irrelevant at this point. 
			for(int particleCounter = 0; particleCounter < (int)newParticles.size(); particleCounter++)
			{
				typedef boost::color_traits<boost::default_color_type> Color;
				std::fill(scratch.colorVector.begin(), scratch.colorVector.end(), Color::white());
				boost::default_dfs_visitor visitor;
				//Check which vertices are accessible from the source or sink, via edges with capacity HIGH_CAPACITY
				boost::detail::depth_first_visit_fixed_impl(undirectedGraph, interestVertices[0], visitor, &(scratch.colorVector[0]), fixedSearchStack, &(newParticles[particleCounter].capacity[0]), boost::detail::nontruth2());
				std::vector<int> newDownEdges;
				for(int i = 0; i < (int)newParticles[particleCounter].downEdges.size(); i++)
				{
					int edgeIndex = newParticles[particleCounter].downEdges[i];
					int vertex1 = edges[edgeIndex].first;
					int vertex2 = edges[edgeIndex].second;
					if(scratch.colorVector[vertex1] == Color::black() && scratch.colorVector[vertex2] == Color::black())
					{
						newParticles[particleCounter].capacity[2 * edgeIndex] = newParticles[particleCounter].capacity[2 * edgeIndex + 1] = HIGH_CAPACITY;
						newParticles[particleCounter].residual[2 * edgeIndex] = newParticles[particleCounter].residual[2 * edgeIndex + 1] = HIGH_CAPACITY;
					}
					else newDownEdges.push_back(edgeIndex);
				}
				std::fill(scratch.colorVector.begin(), scratch.colorVector.end(), Color::white());
				boost::detail::depth_first_visit_fixed_impl(undirectedGraph, interestVertices[1], visitor, &(scratch.colorVector[0]), fixedSearchStack, &(newParticles[particleCounter].capacity[0]), boost::detail::nontruth2());
				for(int i = 0; i < (int)newParticles[particleCounter].downEdges.size(); i++)
				{
					int edgeIndex = newParticles[particleCounter].downEdges[i];
					int vertex1 = edges[edgeIndex].first;
					int vertex2 = edges[edgeIndex].second;
					if(scratch.colorVector[vertex1] == Color::black() && scratch.colorVector[vertex2] == Color::black())
					{
						newParticles[particleCounter].capacity[2 * edgeIndex] = newParticles[particleCounter].capacity[2 * edgeIndex + 1] = HIGH_CAPACITY;
						newParticles[particleCounter].residual[2 * edgeIndex] = newParticles[particleCounter].residual[2 * edgeIndex + 1] = HIGH_CAPACITY;
					}
					else newDownEdges.push_back(edgeIndex);
				}
				newParticles[particleCounter].downEdges.swap(newDownEdges);
			}
			//Sort by state.
			std::sort(newParticles.begin(), newParticles.end(), [nEdges](const approximateZeroVarianceWORMergeImpl::particle& first, const approximateZeroVarianceWORMergeImpl::particle& second){ return memcmp(&(first.capacity[0]), &(second.capacity[0]), sizeof(int)*nEdges*2) < 0;});
			std::size_t unitsAfterMerge = 0;
			//Merge particles
			{
				int particleCounter = 0;
				while(particleCounter < (int)newParticles.size())
				{
					mpfr_class additionalWeight = 0, additionalImportanceDensity = 0;
					int particleCounter2 = particleCounter+1;
					for(; particleCounter2 < (int)newParticles.size() && memcmp(&(newParticles[particleCounter].capacity[0]), &(newParticles[particleCounter2].capacity[0]), sizeof(int)*2 * nEdges) == 0; particleCounter2++)
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
					if(newParticles[i].weight != 0) particles.emplace_back(std::move(newParticles[i]));
				}
			}
			else
			{
				indices.clear();
				sampling::sampfordFromParetoNaive(samplingArgs, args.randomSource);
				for(std::vector<int>::iterator i = indices.begin(); i != indices.end(); i++)
				{
					particles.emplace_back(std::move(newParticles[*i]));
					particles.back().importanceDensity /= rescaledWeights[*i];
					particles.back().weight /= rescaledWeights[*i];
				}
			}
		}
		for(std::vector<approximateZeroVarianceWORMergeImpl::particle>::iterator i = particles.begin(); i != particles.end(); i++)
		{
			args.estimate += i->weight;
		}
	}
}
