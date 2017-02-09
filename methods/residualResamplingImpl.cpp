#include "residualResamplingImpl.h"
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include <boost/random/random_number_generator.hpp>
#include <boost/random/uniform_01.hpp>
#include "aliasMethod.h"
namespace networkReliability
{
	namespace residualResamplingImpl
	{
		struct residualResamplingScratch
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
			std::vector<int> capacity, residual;
			mpfr_class weight;
			int minCutSize;
			particle()
			{}
			particle(const particle& other)
				: capacity(other.capacity), residual(other.residual), weight(other.weight), minCutSize(other.minCutSize)
			{}
			particle& operator=(const particle& other)
			{
				capacity = other.capacity;
				residual = other.residual;
				weight = other.weight;
				minCutSize = other.minCutSize;
				return *this;
			}
			particle(particle&& other)
				: capacity(std::move(other.capacity)), residual(std::move(other.residual)), weight(std::move(other.weight)), minCutSize(other.minCutSize)
			{}
			particle& operator=(particle&& other)
			{
				capacity.swap(other.capacity);
				residual.swap(other.residual);
				weight = other.weight;
				minCutSize = other.minCutSize;
				return *this;
			}
		};
	}
	int getMinCut(std::vector<int>& capacity, std::vector<int>& residual, const context::internalDirectedGraph& graph, const context::internalGraph& undirectedGraph, const std::vector<int>& interestVertices, residualResamplingImpl::residualResamplingScratch& scratch)
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
	void residualResampling(residualResamplingArgs& args)
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

		//Cache powers of the inopProbability 
		boost::scoped_array<mpfr_class> cachedInopPowers(new mpfr_class[nEdges+1]);
		for(std::size_t i = 0; i < nEdges+1; i++)
		{
			cachedInopPowers[i] = boost::multiprecision::pow(inopProbability, i);
		}
		//Temporaries for calculating max flow values
		residualResamplingImpl::residualResamplingScratch scratch;
		scratch.colorVector.resize(nVertices);

		std::vector<std::ptrdiff_t> aliasMethodTemporary1, aliasMethodTemporary2;
		std::vector<std::pair<double, std::ptrdiff_t> > aliasMethodTemporary3;
		std::vector<double> remainders(args.n);

		boost::random::uniform_01<float,float> uniformReal;

		//Initialise with the two initial choices - The first edge can be up or down. 
		std::vector<residualResamplingImpl::particle> particles, newParticles;
		{
			int minCutSizeDown, minCutSizeUp;
			std::vector<int> capacity(2*nEdges, 1), residual(2*nEdges, 1);

			capacity[0] = capacity[1] = 0;
			residual[0] = residual[1] = 0;
			minCutSizeDown = getMinCut(capacity, residual, graph, undirectedGraph, interestVertices, scratch);
			capacity[0] = capacity[1] = HIGH_CAPACITY;
			residual[0] = residual[1] = HIGH_CAPACITY;
			minCutSizeUp = getMinCut(capacity, residual, graph, undirectedGraph, interestVertices, scratch);

			mpfr_class minCutDownProb = cachedInopPowers[minCutSizeDown];
			mpfr_class minCutUpProb = cachedInopPowers[minCutSizeUp];
			mpfr_class qTilde = inopProbability * minCutDownProb;
			qTilde = qTilde / (qTilde  + opProbability * minCutUpProb);
			for(int i = 0; i < (int)args.n; i++)
			{
				residualResamplingImpl::particle newParticle;
				newParticle.capacity.resize(2*nEdges, 1);
				newParticle.residual.resize(2*nEdges, 1);
				if(minCutSizeUp >= HIGH_CAPACITY || uniformReal(args.randomSource) < qTilde.convert_to<double>())
				{
					newParticle.capacity[0] = newParticle.capacity[1] = 0;
					newParticle.residual[0] = newParticle.residual[1] = 0;
					newParticle.weight = inopProbability / qTilde;
					newParticle.minCutSize = minCutSizeDown;
				}
				else
				{
					newParticle.capacity[0] = newParticle.capacity[1] = HIGH_CAPACITY;
					newParticle.residual[0] = newParticle.residual[1] = HIGH_CAPACITY;
					newParticle.weight = opProbability / (1 - qTilde);
					newParticle.minCutSize = minCutSizeUp;
				}
				newParticles.emplace_back(std::move(newParticle));
			}
			particles.swap(newParticles);
		}
		for(int edgeCounter = 1; edgeCounter < (int)nEdges; edgeCounter++)
		{
			for(int particleCounter = 0; particleCounter < (int)particles.size(); particleCounter++)
			{
				residualResamplingImpl::particle& currentParticle = particles[particleCounter];

				currentParticle.capacity[2*edgeCounter] = currentParticle.capacity[2*edgeCounter+1] = 0;
				currentParticle.residual[2*edgeCounter] = currentParticle.residual[2*edgeCounter+1] = 0;
				//int minCutSizeDown = getMinCut(currentParticle.capacity, currentParticle.residual, graph, undirectedGraph, interestVertices, scratch);
				int minCutSizeDown = args.contextObj.getMinCut(currentParticle.capacity);
				currentParticle.capacity[2*edgeCounter] = currentParticle.capacity[2*edgeCounter+1] = HIGH_CAPACITY;
				currentParticle.residual[2*edgeCounter] = currentParticle.residual[2*edgeCounter+1] = HIGH_CAPACITY;
				//int minCutSizeUp = getMinCut(currentParticle.capacity, currentParticle.residual, graph, undirectedGraph, interestVertices, scratch);
				int minCutSizeUp = args.contextObj.getMinCut(currentParticle.capacity);

				mpfr_class minCutDownProb = cachedInopPowers[minCutSizeDown];
				mpfr_class minCutUpProb;
				if(minCutSizeUp >= HIGH_CAPACITY) minCutUpProb = 0;
				else minCutUpProb = cachedInopPowers[minCutSizeUp];
				mpfr_class qTilde = inopProbability * minCutDownProb;
				qTilde = qTilde / (qTilde  + opProbability * minCutUpProb);

				if(minCutSizeUp >= HIGH_CAPACITY)
				{
					currentParticle.capacity[2*edgeCounter] = currentParticle.capacity[2*edgeCounter+1] = 0;
					currentParticle.residual[2*edgeCounter] = currentParticle.residual[2*edgeCounter+1] = 0;
					currentParticle.weight *= inopProbability;
					currentParticle.minCutSize = minCutSizeDown;
				}
				else if(uniformReal(args.randomSource) < qTilde.convert_to<double>())
				{
					currentParticle.capacity[2*edgeCounter] = currentParticle.capacity[2*edgeCounter+1] = 0;
					currentParticle.residual[2*edgeCounter] = currentParticle.residual[2*edgeCounter+1] = 0;
					currentParticle.weight *= inopProbability / qTilde;
					currentParticle.minCutSize = minCutSizeDown;
				}
				else
				{
					currentParticle.capacity[2*edgeCounter] = currentParticle.capacity[2*edgeCounter+1] = HIGH_CAPACITY;
					currentParticle.residual[2*edgeCounter] = currentParticle.residual[2*edgeCounter+1] = HIGH_CAPACITY;
					currentParticle.weight *= opProbability / (1 - qTilde);
					currentParticle.minCutSize = minCutSizeUp;
				}
			}
			newParticles.clear();

			mpfr_class average = 0;
			for(int particleCounter = 0; particleCounter < (int)particles.size(); particleCounter++) average += particles[particleCounter].weight;
			average /= particles.size();

			int remainingParticles = args.n;
			for(int particleCounter = 0; particleCounter < (int)particles.size(); particleCounter++)
			{
				int copies = (int)mpfr_class(particles[particleCounter].weight / average).convert_to<double>();
				for(int i = 0; i < copies; i++) 
				{
					newParticles.push_back(particles[particleCounter]);
					newParticles.back().weight = average;
				}
				remainders[particleCounter] = mpfr_class(particles[particleCounter].weight - copies * average).convert_to<double>();
				remainingParticles -= copies;
			}
			
			mpfr_class remaindersAverage = 0;
			for(int particleCounter = 0; particleCounter < (int)particles.size(); particleCounter++) remaindersAverage += remainders[particleCounter];
			double sumRemainders = remaindersAverage.convert_to<double>();
			remaindersAverage /= remainingParticles;

			aliasMethod::aliasMethod alias(remainders, sumRemainders, aliasMethodTemporary1, aliasMethodTemporary2, aliasMethodTemporary3);
			for (std::size_t k = 0; k < (std::size_t)remainingParticles; k++)
			{
				int index = (int)alias(args.randomSource);
				newParticles.push_back(particles[index]);
				newParticles.back().weight = remaindersAverage;
			}
			particles.swap(newParticles);
		}
		args.estimate = 0;
		for(std::vector<residualResamplingImpl::particle>::iterator i = particles.begin(); i != particles.end(); i++)
		{
			args.estimate += i->weight;
		}
		args.estimate /= args.n;
	}
}
