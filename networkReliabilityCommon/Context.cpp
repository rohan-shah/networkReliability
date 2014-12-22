#include "Context.h"
#include <boost/graph/graphml.hpp>
#include <fstream>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/scoped_array.hpp>
#include <boost/graph/stoer_wagner_min_cut.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include <boost/serialization/map.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include "allPointsMaxFlow.hpp"
namespace networkReliability
{
	namespace ContextImpl
	{
		template<typename Key, Context::internalGraph::vertices_size_type Ret> class constant_property_map_vertices_size_type : public boost::put_get_helper<Context::internalGraph::vertices_size_type, constant_property_map_vertices_size_type<Key, Ret> > 
		{
		public:
			typedef Key key_type;
			typedef Context::internalGraph::vertices_size_type reference;
			typedef Context::internalGraph::vertices_size_type value_type;

			typedef boost::readable_property_map_tag category;

			constant_property_map_vertices_size_type(){}

			reference operator[](const Key&) const 
			{
				return Ret;
			}
		};
		template<typename Key, int Ret> class constant_property_map_int : public boost::put_get_helper<int, constant_property_map_int<Key, Ret> > 
		{
		public:
			typedef Key key_type;
			typedef int reference;
			typedef int value_type;

			typedef boost::readable_property_map_tag category;

			constant_property_map_int(){}

			reference operator[](const Key&) const 
			{
				return Ret;
			}
		};

		struct twoDArray
		{
			int* base;
			std::size_t dim;
			struct twoDArrayInternal
			{
				twoDArrayInternal(int* base)
					:base(base)
				{};
				int& operator[](std::size_t j)
				{
					return *(base + j);
				}
				const int& operator[](std::size_t j) const
				{
					return *(base + j);
				}
				int* base;
			};
			twoDArrayInternal operator[](std::size_t i) const
			{
				return twoDArrayInternal(base + dim*i);
			};
		};
	}
	Context& Context::operator=(Context&& other)
	{
		graph.swap(other.graph);
		interestVertices.swap(other.interestVertices);
		vertexPositions.swap(other.vertexPositions);
		edgeDistances.swap(other.edgeDistances);
		std::swap(operationalProbability, other.operationalProbability);
		allOpDistributions.swap(other.allOpDistributions);
		allInopDistributions.swap(other.allInopDistributions);
		std::swap(minCutEdges, other.minCutEdges);

		edgeResidualCapacityVector.swap(other.edgeResidualCapacityVector);
		capacityVector.swap(other.capacityVector);
		directedGraph.swap(other.directedGraph);
		vertexPredecessorVector.swap(other.vertexPredecessorVector);
		colorVector.swap(other.colorVector);
		distanceVector.swap(other.distanceVector);

		std::swap(nEdges, other.nEdges);
		std::swap(inoperationalProbabilityD, other.inoperationalProbabilityD);
		std::swap(_useMinCut, other._useMinCut);
		return *this;
	}
	Context::Context(Context&& other)
	{
		graph.swap(other.graph);
		interestVertices.swap(other.interestVertices);
		vertexPositions.swap(other.vertexPositions);
		edgeDistances.swap(other.edgeDistances);
		std::swap(operationalProbability, other.operationalProbability);
		allOpDistributions.swap(other.allOpDistributions);
		allInopDistributions.swap(other.allInopDistributions);
		std::swap(minCutEdges, other.minCutEdges);

		edgeResidualCapacityVector.swap(other.edgeResidualCapacityVector);
		capacityVector.swap(other.capacityVector);
		directedGraph.swap(other.directedGraph);
		vertexPredecessorVector.swap(other.vertexPredecessorVector);
		colorVector.swap(other.colorVector);
		distanceVector.swap(other.distanceVector);

		std::swap(nEdges, other.nEdges);
		std::swap(inoperationalProbabilityD, other.inoperationalProbabilityD);
		std::swap(_useMinCut, other._useMinCut);
	}
	Context::Context(boost::shared_ptr<const inputGraph> unorderedGraph, boost::shared_ptr<const std::vector<int> > edgeOrdering, boost::shared_ptr<const std::vector<int> > interestVertices, boost::shared_ptr<std::vector<vertexPosition> > vertexPositions, const mpfr_class& operationalProbability, boost::shared_array<double> inputEdgeDistances)
		:vertexPositions(vertexPositions), interestVertices(interestVertices), operationalProbability(operationalProbability), _useMinCut(false)
	{
		mpfr_class inoperationalProbability = (1 - operationalProbability);
		inoperationalProbabilityD = inoperationalProbability.convert_to<double>();

		std::size_t nVertices = boost::num_vertices(*unorderedGraph);
		nEdges = boost::num_edges(*unorderedGraph);
		
		edgeResidualCapacityVector.resize(2*nEdges);
		vertexPredecessorVector.resize(2*nEdges);
		capacityVector.resize(2*nEdges, 1);
		vertexPredecessorVector.resize(2*nEdges);
		colorVector.resize(2*nEdges);
		distanceVector.resize(2*nEdges);

		if(nEdges != edgeOrdering->size())
		{
			throw std::runtime_error("Graph ordering data had the wrong size");
		}

		if(nVertices != vertexPositions->size())
		{
			throw std::runtime_error("Vertex position data had the wrong size");
		}
		if(*std::min_element(edgeOrdering->begin(), edgeOrdering->end()) != 0)
		{
			throw std::runtime_error("Wrong minimum vertex in ordering");
		}
		if(*std::max_element(edgeOrdering->begin(), edgeOrdering->end()) != nEdges-1)
		{
			throw std::runtime_error("Wrong maximum vertex in ordering");
		}

		boost::shared_ptr<internalGraph> orderedGraph(new internalGraph(nVertices));
		inputGraph::edge_iterator start, end;
		boost::tie(start, end) = boost::edges(*unorderedGraph);
		int index = 0;
		for(; start != end; start++)
		{
			boost::add_edge(start->m_source, start->m_target, (*edgeOrdering)[index], *orderedGraph);
			index++;
		}

		graph = orderedGraph;
		constructDirectedGraph();
		if(inputEdgeDistances)
		{
			edgeDistances = inputEdgeDistances;
		}
		else
		{
			constructEdgeDistances();
		}

		//are we looking at the all-terminal reliability problem?
		if(interestVertices->size() == nVertices)
		{
			ContextImpl::constant_property_map_vertices_size_type<Context::internalGraph::edge_descriptor, 1L> edgeWeights;

			//BOOST_AUTO(parities, boost::make_one_bit_color_map(num_vertices(*graph), get(boost::vertex_index, *graph)));
			minCutEdges = boost::stoer_wagner_min_cut(*graph, edgeWeights);//, boost::parity_map(parities));
		}
		//or are we looking at the 2-terminal reliability problem?
		else if(interestVertices->size() == 2)
		{
			typedef boost::property_map<Context::internalDirectedGraph, boost::edge_index_t>::const_type edgeIndexMap;
			typedef boost::iterator_property_map<typename std::vector<int>::iterator, edgeIndexMap> edgeCapacityMap;

			edgeIndexMap edgeIndices = boost::get(boost::edge_index, *directedGraph);
			edgeCapacityMap residualCapacityMap(edgeResidualCapacityVector.begin(), edgeIndices);
			edgeCapacityMap capacityMap(capacityVector.begin(), edgeIndices);

			minCutEdges = boost::push_relabel_max_flow(*directedGraph, (*interestVertices)[0], (*interestVertices)[1], boost::residual_capacity_map(residualCapacityMap).capacity_map(capacityMap));
		}
		else
		{
			std::vector<int> flowMatrix(nVertices*nVertices);
			allPointsMaxFlow::allPointsMaxFlowScratch<Context::internalDirectedGraph, int> scratch;
			allPointsMaxFlow::allPointsMaxFlow<Context::internalDirectedGraph, int>(flowMatrix, capacityVector, *directedGraph.get(), scratch);
			int maxFlow = std::numeric_limits<int>::max();
			for(int i = 0; i < interestVertices->size(); i++)
			{
				for(int j = i+1; j < interestVertices->size(); j++)
				{
					maxFlow = std::min(maxFlow, flowMatrix[(*interestVertices)[i] + (*interestVertices)[j] * nVertices]);
				}
			}
			minCutEdges = maxFlow;
		}
	}
	void Context::constructEdgeDistances()
	{
		const std::size_t nEdges = boost::num_edges(*graph);
		const std::size_t nVertices = boost::num_vertices(*graph);
		boost::scoped_array<int> vertexDistances(new int[nVertices * nVertices]);
		int* vertexDistancePtr = vertexDistances.get();

		ContextImpl::twoDArray tmp;
		tmp.base = vertexDistances.get();
		tmp.dim = nVertices;

		ContextImpl::constant_property_map_int<Context::inputGraph::edge_descriptor, 1> edgeWeights;
		boost::johnson_all_pairs_shortest_paths(*graph, tmp, boost::weight_map(edgeWeights));
		
		
		edgeDistances = boost::shared_array<double>(new double[nEdges * nEdges]);
		double* edgeDistancePtr = edgeDistances.get();
		memset(edgeDistancePtr, (int)(sizeof(int)*nEdges*nEdges), 0);

		Context::internalGraph::edge_iterator currentFirst, end;
		boost::tie(currentFirst, end) = boost::edges(*graph);
		while(currentFirst != end)
		{
			int firstIndex = boost::get(boost::edge_index, *graph, *currentFirst);

			Context::internalGraph::edge_iterator currentSecond = currentFirst;
			currentSecond++;
			while(currentSecond != end)
			{
				int secondIndex = boost::get(boost::edge_index, *graph, *currentSecond);
				int possibilities[4] = 
				{
					vertexDistancePtr[currentFirst->m_source + nVertices*currentSecond->m_source], 
					vertexDistancePtr[currentFirst->m_source + nVertices*currentSecond->m_target], 
					vertexDistancePtr[currentFirst->m_target + nVertices*currentSecond->m_target], 
					vertexDistancePtr[currentFirst->m_target + nVertices*currentSecond->m_source]
				};
				
				edgeDistancePtr[secondIndex + nEdges * firstIndex] = edgeDistancePtr[firstIndex + nEdges * secondIndex] = *std::min_element(possibilities + 0, possibilities + 4)+1;
				currentSecond++;
			}
			currentFirst++;
		}
	}
	Context::Context()
		:graph(), directedGraph(), vertexPositions(), nEdges(0), _useMinCut(false)
	{}
	Context Context::emptyContext()
	{
		Context result;

		boost::shared_ptr<Context::internalGraph> graph(new Context::internalGraph(2));
		boost::add_edge(0, 1, 0, *graph);
		result.graph = boost::static_pointer_cast<const Context::internalGraph>(graph);

		boost::shared_ptr<Context::internalDirectedGraph> directedGraph(new Context::internalDirectedGraph(2));
		boost::add_edge(0, 1, 0, *directedGraph);
		boost::add_edge(1, 0, 1, *directedGraph);
		result.directedGraph = boost::static_pointer_cast<const Context::internalDirectedGraph>(directedGraph);

		boost::shared_ptr<std::vector<int> > interestVertices(new std::vector<int>(2));
		(*interestVertices)[0] = 0; (*interestVertices)[1] = 1;
		result.interestVertices = boost::static_pointer_cast<const std::vector<int> >(interestVertices);

		boost::shared_ptr<std::vector<vertexPosition> > vertexPositions(new std::vector<vertexPosition>(2));
		(*vertexPositions)[0] = vertexPosition(vertexPosition::first_type(0.0), vertexPosition::second_type(0.0)); (*vertexPositions)[1] = vertexPosition(vertexPosition::first_type(10.0), vertexPosition::second_type(0.0));
		result.vertexPositions = boost::static_pointer_cast<const std::vector<vertexPosition> >(vertexPositions);
		
		boost::shared_array<double> edgeDistances(new double[4]);
		edgeDistances[0] = edgeDistances[3] = 0;
		edgeDistances[1] = edgeDistances[2] = 1;
		result.edgeDistances.swap(edgeDistances);

		result.operationalProbability = 0.5;
		result.inoperationalProbabilityD = 0.5;
		
		result.minCutEdges = 1;
		result.edgeResidualCapacityVector.resize(2);
		result.capacityVector.resize(2);
		result.vertexPredecessorVector.resize(2);
		result.colorVector.resize(2);
		result.distanceVector.resize(2);

		result.nEdges = 0;
		result._useMinCut = false;

		return result;
	}
	Context Context::gridContext(int gridDimension, boost::shared_ptr<const std::vector<int> > interestVertices, const mpfr_class& operationalProbability)
	{
		boost::shared_ptr<std::vector<vertexPosition> > vertexPositions(new std::vector<vertexPosition>(gridDimension * gridDimension));
		boost::shared_ptr<Context::inputGraph> graph(new Context::inputGraph(gridDimension * gridDimension));
		for(int i = 0; i < gridDimension; i++)
		{
			for(int j = 0; j < gridDimension; j++)
			{
				(*vertexPositions)[i + j * gridDimension] = vertexPosition((float)i*100, (float)j*100);
				if(i != gridDimension - 1) boost::add_edge(i + j*gridDimension, i + 1 +j*gridDimension, *graph);
				if(j != gridDimension - 1) boost::add_edge(i + j*gridDimension, i + (j+1)*gridDimension, *graph);
			}
		}
		boost::shared_ptr<std::vector<int> > edgeOrdering(new std::vector<int>(boost::num_edges(*graph)));
		for(int i = 0; i < edgeOrdering->size(); i++) (*edgeOrdering)[i] = i;
		return Context(graph, edgeOrdering, interestVertices, vertexPositions, operationalProbability);
	}
	Context Context::completeContext(int nVertices, int nInterestVertices, const mpfr_class& operationalProbability)
	{
		if(nInterestVertices < 2 || nInterestVertices > nVertices)
		{
			throw std::runtime_error("Internal error");
		}
		boost::shared_ptr<std::vector<vertexPosition> > vertexPositions(new std::vector<vertexPosition>(nVertices));
		boost::shared_ptr<Context::inputGraph> graph(new Context::inputGraph(nVertices));

		boost::shared_ptr<std::vector<int> > interestVertices(new std::vector<int>(nInterestVertices));
		std::copy(boost::counting_iterator<int>(1), boost::counting_iterator<int>(nInterestVertices), (*interestVertices).begin());

		const double pi = 3.14159265359;
		for(int i = 0; i < nVertices; i++)
		{
			(*vertexPositions)[i] = vertexPosition((float)cos(2*pi*i/nVertices), (float)sin(2*pi*i/nVertices));
			for(int j = i+1; j < nVertices; j++)
			{
				boost::add_edge(i, j, *graph);
			}
		}
		boost::shared_ptr<std::vector<int> > edgeOrdering(new std::vector<int>(boost::num_edges(*graph)));
		for(int i = 0; i < edgeOrdering->size(); i++) (*edgeOrdering)[i] = i;

		return Context(graph, edgeOrdering, interestVertices, vertexPositions, operationalProbability);
	}
	const Context::internalGraph& Context::getGraph() const
	{
		return *graph;
	}
	const std::vector<int>& Context::getInterestVertices() const
	{
		return *interestVertices;
	}
	const std::vector<Context::vertexPosition>& Context::getVertexPositions() const
	{
		return *vertexPositions;
	}
	Context Context::fromFile(std::string path, bool& successful, boost::shared_ptr<const std::vector<int> > interestVertices, std::string& message, const mpfr_class& operationalProbability, bool useSpatialDistances)
	{
		std::ifstream input(path.c_str());
		if(!input.is_open())
		{
			successful = false;
			return Context();
		}
		boost::dynamic_properties properties;
		boost::shared_ptr<inputGraph> graph(new inputGraph());

		typedef boost::property_map<Context::inputGraph, boost::edge_index_t>::type edgeIndexMapType; 
		edgeIndexMapType orderingProperty;
		properties.property("order", orderingProperty);

		boost::vector_property_map<float> xProperty, yProperty;
		properties.property("x", xProperty);
		properties.property("y", yProperty);

		try
		{
			boost::read_graphml(input, *graph, properties);
		}
		catch(...)
		{
			message = "Error parsing graphml file";
			successful = false;
			return Context();
		}
		const std::size_t nEdges = boost::num_edges(*graph);

		boost::shared_ptr<std::vector<vertexPosition> > vertexPositions(new std::vector<vertexPosition>());
		for(auto xIterator = xProperty.storage_begin(), yIterator = yProperty.storage_begin(); xIterator != xProperty.storage_end(); xIterator++, yIterator++)
		{
			vertexPositions->push_back(vertexPosition(*xIterator, *yIterator));
		}

		successful = true;
		//check that the edge_index values are either all 0, or distinct. 
		//In the process copy the values into a vector
		bool allZero = true;
		boost::shared_ptr<std::vector<int> > ordering(new std::vector<int>(boost::num_edges(*graph)));
		Context::inputGraph::edge_iterator current, end;
		boost::tie(current, end) = boost::edges(*graph);
		for (; current != end; current++)
		{
			if (boost::get(boost::edge_index, *graph, *current) != 0)
			{
				allZero = false;
				break;
			}
		}
		if (allZero)
		{
			std::copy(boost::counting_iterator<int>(0), boost::counting_iterator<int>((int)nEdges), ordering->begin());
		}
		else
		{
			std::vector<bool> foundEdgeIndex(nEdges, false);
			boost::tie(current, end) = boost::edges(*graph);
			std::vector<int>::iterator outputIterator = ordering->begin();
			for (; current != end; current++)
			{
				int edgeIndex = boost::get(boost::edge_index, *graph, *current);
				if (edgeIndex >= nEdges || edgeIndex < 0)
				{
					message = "Invalid order value specified";
					successful = false;
					return Context();
				}
				if (foundEdgeIndex[edgeIndex])
				{
					successful = false;
					message = "Duplicate values in ordering";
					return Context();
				}
				else foundEdgeIndex[edgeIndex] = true;
				*outputIterator = edgeIndex;
				outputIterator++;
			}
		}


		int maxInterest = *std::max_element(interestVertices->begin(), interestVertices->end());
		int minInterest = *std::min_element(interestVertices->begin(), interestVertices->end());
		if(minInterest < 0 || maxInterest >= boost::num_vertices(*graph))
		{
			successful = false;
			message = "Invalid vertex indices entered for input interestVertices";
			return Context();
		}
		boost::shared_array<double> inputEdgeDistances;
		if(useSpatialDistances)
		{
			inputEdgeDistances.reset(new double[nEdges*nEdges]);
			memset(inputEdgeDistances.get(), 0, sizeof(double)*nEdges*nEdges);
			Context::inputGraph::edge_iterator end, firstEdge, secondEdge;

			boost::tie(firstEdge, end) = boost::edges(*graph);
			std::vector<vertexPosition> edgeCentres(nEdges);
			std::vector<int>::iterator orderIterator = ordering->begin();
			for(;firstEdge != end; firstEdge++,orderIterator++)
			{
				int edgeIndex = *orderIterator;
				vertexPosition& currentCentre = edgeCentres[edgeIndex];
				currentCentre.first = ((*vertexPositions)[firstEdge->m_source].first + (*vertexPositions)[firstEdge->m_target].first)/2;
				currentCentre.second = ((*vertexPositions)[firstEdge->m_source].second + (*vertexPositions)[firstEdge->m_target].second)/2;
			}

			boost::tie(firstEdge, end) = boost::edges(*graph);
			orderIterator = ordering->begin();
			for(;firstEdge != end; firstEdge++, orderIterator++)
			{
				int firstEdgeIndex = *orderIterator;
				boost::tie(secondEdge, end) = boost::edges(*graph);
				vertexPosition& firstEdgeCentre = edgeCentres[firstEdgeIndex];
				std::vector<int>::iterator secondOrderIterator = ordering->begin();
				for(; secondEdge != end; secondEdge++, secondOrderIterator++)
				{
					if(firstEdge != secondEdge)
					{
						int secondEdgeIndex = *secondOrderIterator;
						vertexPosition& secondEdgeCentre = edgeCentres[secondEdgeIndex];
						inputEdgeDistances[firstEdgeIndex + secondEdgeIndex * nEdges] = sqrt((firstEdgeCentre.first - secondEdgeCentre.first)*(firstEdgeCentre.first - secondEdgeCentre.first) + (firstEdgeCentre.second - secondEdgeCentre.second)*(firstEdgeCentre.second - secondEdgeCentre.second));
					}
				}
			}
		}
		return Context(graph, ordering, interestVertices, vertexPositions, operationalProbability, inputEdgeDistances);
	}
	const double* Context::getEdgeDistances() const
	{
		return edgeDistances.get();
	}
	const mpfr_class& Context::getOperationalProbability() const
	{
		return operationalProbability;
	}
	Context::~Context()
	{
	}
	const ::TruncatedBinomialDistribution::TruncatedBinomialDistribution& Context::getOpDistribution(std::size_t firstAllowedValue, std::size_t lastAllowedValue, std::size_t n) const
	{
		::TruncatedBinomialDistribution::TruncatedBinomialDistribution::key key;
		key.firstAllowedValue = firstAllowedValue;
		key.lastAllowedValue = lastAllowedValue;
		key.n = n;

		::TruncatedBinomialDistribution::TruncatedBinomialDistributionCollection::mapType::iterator searchIterator;
		searchIterator = allOpDistributions.data.find(key);
		if (searchIterator == allOpDistributions.data.end())
		{
			::TruncatedBinomialDistribution::TruncatedBinomialDistributionCollection::mapType::value_type toInsert(key, ::TruncatedBinomialDistribution::TruncatedBinomialDistribution(n, firstAllowedValue, lastAllowedValue, operationalProbability));
			return allOpDistributions.data.insert(std::move(toInsert)).first->second;
		}
		else return searchIterator->second;
	}
	const ::TruncatedBinomialDistribution::TruncatedBinomialDistribution& Context::getInopDistribution(std::size_t firstAllowedValue, std::size_t lastAllowedValue, std::size_t n) const
	{
		::TruncatedBinomialDistribution::TruncatedBinomialDistribution::key key;
		key.firstAllowedValue = firstAllowedValue;
		key.lastAllowedValue = lastAllowedValue;
		key.n = n;

		::TruncatedBinomialDistribution::TruncatedBinomialDistributionCollection::mapType::iterator searchIterator;
		searchIterator = allInopDistributions.data.find(key);
		if(searchIterator == allInopDistributions.data.end())
		{
			::TruncatedBinomialDistribution::TruncatedBinomialDistributionCollection::mapType::value_type toInsert(key, ::TruncatedBinomialDistribution::TruncatedBinomialDistribution(n, firstAllowedValue, lastAllowedValue, 1-operationalProbability));
			return allInopDistributions.data.insert(std::move(toInsert)).first->second;
		}
		else return searchIterator->second;
	}
	std::size_t Context::getMinCutEdges() const
	{
		return minCutEdges;
	}
	void Context::constructDirectedGraph()
	{
		const std::size_t nEdges = boost::num_edges(*graph);
		boost::shared_ptr<internalDirectedGraph> directedGraph(new internalDirectedGraph(boost::num_vertices(*graph)));
		internalDirectedGraph& directedGraphRef = *directedGraph;
	
		boost::property_map<internalDirectedGraph, boost::edge_reverse_t>::type reverseMap = boost::get(boost::edge_reverse, directedGraphRef);

		internalGraph::edge_iterator start, end, current;
		boost::tie(start, end) = boost::edges(*graph);
		//add edges
		current = start;
		while(current != end)
		{
			int originalEdgeIndex = boost::get(boost::edge_index, *graph, *current);
			std::pair<internalDirectedGraph::edge_descriptor, bool> firstEdgePair = boost::add_edge(current->m_source, current->m_target, 2*originalEdgeIndex, directedGraphRef);
			std::pair<internalDirectedGraph::edge_descriptor, bool> secondEdgePair = boost::add_edge(current->m_target, current->m_source, 2*originalEdgeIndex+1, directedGraphRef);

			boost::put(reverseMap, firstEdgePair.first, secondEdgePair.first);
			boost::put(reverseMap, secondEdgePair.first, firstEdgePair.first);

			current++;
		}

		this->directedGraph = boost::static_pointer_cast<const internalDirectedGraph>(directedGraph);
	}
	const Context::internalDirectedGraph& Context::getDirectedGraph() const
	{
		return *directedGraph;
	}
	std::vector<int>& Context::getCapacityVector() const
	{
		return capacityVector;
	}
	int Context::getMinCut(std::vector<int>& capacityVector) const
	{
		std::size_t nVertices = boost::num_vertices(*graph);
		//Are we looking at the all-terminal reliability problem?
		if(interestVertices->size() == nVertices)
		{
			ContextImpl::constant_property_map_vertices_size_type<Context::internalGraph::edge_descriptor, 1L> edgeWeights;

			//BOOST_AUTO(parities, boost::make_one_bit_color_map(num_vertices(*graph), get(boost::vertex_index, *graph)));
			return (int)boost::stoer_wagner_min_cut(*graph, edgeWeights);//, boost::parity_map(parities));
		}
		//or are we looking at the k-terminal reliability problem?
		else
		{
			const std::size_t nInterestVertices = interestVertices->size();
			//Use the all-points max flow
			if(interestVertices->size() * (interestVertices->size() - 1) / 2 > nVertices - 1)
			{
				std::vector<int> maxFlowResults(nVertices * nVertices, std::numeric_limits<int>::max());
				allPointsMaxFlow::allPointsMaxFlowScratch<Context::internalDirectedGraph, int> scratch;
				allPointsMaxFlow::allPointsMaxFlow<Context::internalDirectedGraph, int>(maxFlowResults, capacityVector, *directedGraph, scratch);
				int minimum = std::numeric_limits<int>::max();
				for(std::size_t i = 0; i < nInterestVertices; i++)
				{
					for(std::size_t j = i+1; j < nInterestVertices; j++)
					{
						minimum = std::min(minimum, maxFlowResults[(*interestVertices)[i] + (*interestVertices)[j] * nVertices]);
					}
				}
				return minimum;
			}
			//Otherwise just call max-flow a bunch of times (the naive version)
			else
			{
				std::vector<int> maxFlowResults(nInterestVertices * nInterestVertices, std::numeric_limits<int>::max());
				typedef boost::property_map<Context::internalDirectedGraph, boost::edge_index_t>::const_type edgeIndexMapType;
				typedef boost::property_map<Context::internalDirectedGraph, boost::vertex_index_t>::const_type vertexIndexMapType;
				typedef boost::iterator_property_map<typename std::vector<int>::iterator, edgeIndexMapType> edgeCapacityMapType;
				typedef boost::iterator_property_map<typename std::vector<Context::internalDirectedGraph::edge_descriptor>::iterator, vertexIndexMapType> vertexPredecessorMapType;
				typedef boost::iterator_property_map<typename std::vector<boost::default_color_type>::iterator, vertexIndexMapType> colorMapType;
				typedef boost::iterator_property_map<typename std::vector<int>::iterator, vertexIndexMapType> distanceMapType;

				edgeIndexMapType edgeIndexMap = boost::get(boost::edge_index, *directedGraph);
				vertexIndexMapType vertexIndexMap = boost::get(boost::vertex_index, *directedGraph);
				edgeCapacityMapType residualCapacityMap(edgeResidualCapacityVector.begin(), edgeIndexMap);
				edgeCapacityMapType edgeCapacityMap(capacityVector.begin(), edgeIndexMap);
				vertexPredecessorMapType vertexPredecessorMap(vertexPredecessorVector.begin(), vertexIndexMap);
				colorMapType colorMap(colorVector.begin(), vertexIndexMap);
				distanceMapType distanceMap(distanceVector.begin(), vertexIndexMap);

				for(std::size_t i = 0; i < nInterestVertices; i++)
				{
					for(std::size_t j = i+1; j < nInterestVertices; j++)
					{
						maxFlowResults[i + j * nInterestVertices] = maxFlowResults[j + i * nInterestVertices] = boost::boykov_kolmogorov_max_flow(*directedGraph, (*interestVertices)[i], (*interestVertices)[j], boost::residual_capacity_map(residualCapacityMap).capacity_map(edgeCapacityMap).predecessor_map(vertexPredecessorMap).color_map(colorMap).distance_map(distanceMap));
					}
				}
				return *std::min_element(maxFlowResults.begin(), maxFlowResults.end());
			}
		}
	}
	bool Context::useMinCut() const
	{
		return _useMinCut;
	}
	std::size_t Context::getNEdges() const
	{
		return nEdges;
	}
	double Context::getInoperationalProbabilityD() const
	{
		return inoperationalProbabilityD;
	}
	void Context::setMinCut(bool useMinCut)
	{
		this->_useMinCut = useMinCut;
	}
}
