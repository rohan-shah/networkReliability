#ifndef CONTEXT_HEADER_GUARD
#define CONTEXT_HEADER_GUARD
#include <boost/graph/adjacency_list.hpp>
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <vector>
#include "TruncatedBinomialDistribution.h"
#define HIGH_CAPACITY 100000
namespace networkReliability
{
	class Context : public boost::noncopyable
	{
	public:
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

		Context(boost::shared_ptr<const inputGraph> graph, boost::shared_ptr<const std::vector<int> > edgeOrdering, boost::shared_ptr<const std::vector<int> > interestVertices, boost::shared_ptr<std::vector<vertexPosition> > vertexPositions, const mpfr_class& operationalProbability);

		Context& operator=(Context&& other);
		Context(Context&& other);
		~Context();
		
		static Context gridContext(int gridDimension, boost::shared_ptr<const std::vector<int> > interestVertices, const mpfr_class& operationalProbability);
		static Context fromFile(std::string path, bool& successful, boost::shared_ptr<const std::vector<int> > interestVertices, std::string& message, const mpfr_class& operationalProbability);
		static Context emptyContext();
		const internalGraph& getGraph() const;
		const internalDirectedGraph& getDirectedGraph() const;
		const int* getEdgeDistances() const;
		std::size_t getNEdges() const;
		const std::vector<int>& getInterestVertices() const;
		const std::vector<vertexPosition>& getVertexPositions() const;
		const mpfr_class& getOperationalProbability() const;
		void loadDistributions(std::string path);
		const ::TruncatedBinomialDistribution::TruncatedBinomialDistribution& getDistribution(std::size_t firstAllowedValue, std::size_t lastAllowedValue, std::size_t n) const;
		std::size_t getMinCutEdges() const;
		std::vector<int>& getCapacityVector() const;
		int getMinCut(std::vector<int>& capacityVector) const;
		static Context completeContext(int nVertices, int nInterestVertices, const mpfr_class& operationalProbability);
		double getInoperationalProbabilityD() const;
	private:
		//Context& operator=(Context const& other);
		Context();
		void constructEdgeDistances();
		void constructDirectedGraph();
		boost::shared_ptr<const internalGraph> graph;
		boost::shared_ptr<const internalDirectedGraph> directedGraph;
		boost::shared_ptr<const std::vector<int> > interestVertices;
		boost::shared_ptr<const std::vector<vertexPosition> > vertexPositions;

		std::size_t nEdges;

		boost::shared_array<int> edgeDistances;
		mpfr_class operationalProbability;
		double inoperationalProbabilityD;
		mutable ::TruncatedBinomialDistribution::TruncatedBinomialDistributionCollection allDistributions;
		
		std::size_t minCutEdges;

		std::vector<std::string> distributionPaths;
		//These are used in the min paths call. Stored here so they can be reused. 
		mutable std::vector<int> edgeResidualCapacityVector;
		mutable std::vector<int> capacityVector;
		mutable std::vector<internalDirectedGraph::edge_descriptor> vertexPredecessorVector;
		mutable std::vector<boost::default_color_type> colorVector;
		mutable std::vector<int> distanceVector;
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