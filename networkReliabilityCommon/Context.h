#ifndef CONTEXT_HEADER_GUARD
#define CONTEXT_HEADER_GUARD
#include <boost/graph/adjacency_list.hpp>
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/graph/adj_list_serialize.hpp>
#include <vector>
#include "TruncatedBinomialDistribution.h"
#include "allPointsMaxFlow.hpp"
#include "serializeGMP.hpp"
#define HIGH_CAPACITY 100000
namespace networkReliability
{
	class NetworkReliabilitySubObs;
	class Context : private boost::noncopyable
	{
	public:
		Context(boost::archive::binary_iarchive& ar);
		Context(boost::archive::text_iarchive& ar);
		friend class boost::serialization::access;
		typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, boost::property<boost::edge_index_t, int> > inputGraph;
		typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, boost::property<boost::edge_index_t, int> > internalGraph;
		typedef boost::adjacency_list_traits<boost::vecS, boost::vecS, boost::bidirectionalS> directedGraphTraits;
		typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, boost::no_property, 
			boost::property<boost::edge_index_t, int, 
				boost::property<boost::edge_reverse_t, directedGraphTraits::edge_descriptor>
			> 
		> internalDirectedGraph;
		typedef std::pair<float, float> vertexPosition;

		Context(boost::shared_ptr<const inputGraph> graph, boost::shared_ptr<const std::vector<unsigned int> > edgeOrdering, boost::shared_ptr<const std::vector<int> > interestVertices, boost::shared_ptr<std::vector<vertexPosition> > vertexPositions, const mpfr_class& operationalProbability, boost::shared_array<double> inputEdgeDistances = boost::shared_array<double>());

		Context& operator=(Context&& other);
		Context(Context&& other);
		~Context();
		
		static Context gridContext(int gridDimension, boost::shared_ptr<const std::vector<int> > interestVertices, const mpfr_class& operationalProbability);
		static Context fromFile(std::string path, bool& successful, boost::shared_ptr<const std::vector<int> > interestVertices, std::string& message, const mpfr_class& operationalProbability, bool useSpatialDistances);
		static Context emptyContext();
		const internalGraph& getGraph() const;
		const internalDirectedGraph& getDirectedGraph() const;
		const double* getEdgeDistances() const;
		std::size_t getNEdges() const;
		const std::vector<int>& getInterestVertices() const;
		const std::vector<vertexPosition>& getVertexPositions() const;
		const mpfr_class& getOperationalProbability() const;
		const ::TruncatedBinomialDistribution::TruncatedBinomialDistribution& getInopDistribution(std::size_t firstAllowedValue, std::size_t lastAllowedValue, std::size_t n) const;
		const ::TruncatedBinomialDistribution::TruncatedBinomialDistribution& getOpDistribution(std::size_t firstAllowedValue, std::size_t lastAllowedValue, std::size_t n) const;
		std::size_t getMinCutEdges() const;
		std::vector<int>& getCapacityVector() const;
		int getMinCut(std::vector<int>& capacityVector) const;
		static Context completeContext(int nVertices, int nInterestVertices, const mpfr_class& operationalProbability);
		double getInoperationalProbabilityD() const;
		bool useMinCut() const;
		void setMinCut(bool useMinCut);
	private:
		BOOST_SERIALIZATION_SPLIT_MEMBER()
		template<class Archive> void save(Archive& ar, const unsigned int version) const
		{
			std::string typeString = "context";
			ar << typeString;
			//save all shared_ptrs as references instead. 
			ar & *graph;
			///I would have to write extra serialization code to get this line to compile
			//ar & directedGraph;
			ar & *interestVertices;
			ar & *vertexPositions;
			ar & nEdges;
			ar << boost::serialization::make_array(edgeDistances.get(), nEdges*nEdges);
			ar << operationalProbability << inoperationalProbabilityD << minCutEdges;
			typeString = "end_context";
			ar << typeString;
		}
		template<class Archive> void load(Archive& ar, const unsigned int version)
		{
			std::string typeString;
			ar >> typeString;
			if(typeString != "context")
			{
				throw std::runtime_error("Incorrect type specifier");
			}
			//load all shared_ptrs as references
			//have to load as non-const
			{
				boost::shared_ptr<internalGraph> nonConstGraph(new internalGraph());
				ar & *nonConstGraph;
				graph = nonConstGraph;
			}
			//I would have to write extra serialization code to get this line to compile, so just call constructDerictedGraph instead
			//ar & directedGraph;
			constructDirectedGraph();
			//have to load as non-const
			{
				boost::shared_ptr<std::vector<int> > nonConstInterestVertices(new std::vector<int>());
				ar & *nonConstInterestVertices;
				interestVertices = nonConstInterestVertices;
			}
			//have to load as non-const
			{
				boost::shared_ptr<std::vector<vertexPosition> > nonConstVertexPositions(new std::vector<vertexPosition>());
				ar & *nonConstVertexPositions;
				vertexPositions = nonConstVertexPositions;
			}
			ar & nEdges;
			if(boost::num_edges(*graph) != nEdges)
			{
				throw std::runtime_error("Invalid number of edges for graph");
			}
			if(boost::num_edges(*directedGraph) != 2*nEdges)
			{
				throw std::runtime_error("Involid number of edges for directedGraph");
			}
			edgeDistances.reset(new double[nEdges*nEdges]);
			ar >> boost::serialization::make_array(edgeDistances.get(), nEdges*nEdges);
			ar >> operationalProbability >> inoperationalProbabilityD >> minCutEdges;
			ar >> typeString;
			if(typeString != "end_context")
			{
				throw std::runtime_error("File did not end with correct type specifier");
			}
			capacityVector.clear();
			capacityVector.resize(2*nEdges, 1);
			edgeResidualCapacityVector.resize(2*nEdges);
			vertexPredecessorVector.resize(2*nEdges);
			colorVector.resize(2*nEdges);
			distanceVector.resize(2*nEdges);
		}
		bool _useMinCut;
		Context();
		void constructEdgeDistances();
		void constructDirectedGraph();
		boost::shared_ptr<const internalGraph> graph;
		boost::shared_ptr<const internalDirectedGraph> directedGraph;
		boost::shared_ptr<const std::vector<int> > interestVertices;
		boost::shared_ptr<const std::vector<vertexPosition> > vertexPositions;

		std::size_t nEdges;

		boost::shared_array<double> edgeDistances;
		mpfr_class operationalProbability;
		double inoperationalProbabilityD;
		mutable ::TruncatedBinomialDistribution::TruncatedBinomialDistributionCollection allInopDistributions;
		mutable ::TruncatedBinomialDistribution::TruncatedBinomialDistributionCollection allOpDistributions;
		
		std::size_t minCutEdges;

		//These are used in the min paths call. Stored here so they can be reused. 
		mutable std::vector<int> edgeResidualCapacityVector;
		mutable std::vector<int> capacityVector;
		mutable std::vector<internalDirectedGraph::edge_descriptor> vertexPredecessorVector;
		mutable std::vector<boost::default_color_type> colorVector;
		mutable std::vector<int> distanceVector;
		//Temporary storage for Context::getMinCut()
		mutable std::vector<int> maxFlowResults;
		mutable allPointsMaxFlow::allPointsMaxFlowScratch<Context::internalDirectedGraph, int> scratch;
	};
}
namespace boost
{
	namespace networkReliabilityImpl
	{
		struct customStruct
		{};
	}
	template<> struct property_map<const networkReliability::Context::internalDirectedGraph, edge_capacity_t>
	{
		typedef networkReliabilityImpl::customStruct const_type;
	};
	template<> struct property_traits<networkReliabilityImpl::customStruct>
	{
		typedef int value_type;
	};
}
#endif
