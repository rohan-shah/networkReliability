#ifndef NETWORK_RELIABILITY_SUB_OBS_SUB_OBS_ROOT_HEADER_GUARD
#define NETWORK_RELIABILITY_SUB_OBS_SUB_OBS_ROOT_HEADER_GUARD
#include <boost/random/mersenne_twister.hpp>
#include "context.h"
#include "serializeGMP.hpp"
#include "edgeState.h"
#include <boost/shared_array.hpp>
#include "includeMPFR.h"
#include "graphAlgorithms.h"
#include "networkReliabilityObs.h"
namespace networkReliability
{
	namespace subObs
	{
		class subObs : public ::networkReliability::NetworkReliabilityObs
		{
		public:
			void getReducedGraph(Context::internalGraph& outputGraph, std::vector<int>& edgeCounts, std::vector<int>& components, boost::detail::depth_first_visit_restricted_impl_helper<Context::internalGraph>::stackType& stack, std::vector<boost::default_color_type>& colorMap) const;
			typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::property<boost::vertex_name_t, int>, boost::property<boost::edge_index_t, int, boost::property<boost::edge_inop_probability_t, mpfr_class, boost::property<boost::edge_op_probability_t, mpfr_class> > > > reducedGraphWithProbabilities;
			struct getReducedGraphNoSelfWithWeightsInput
			{
			public:
				getReducedGraphNoSelfWithWeightsInput(const std::vector<int>& interestVertices)
					:interestVertices(interestVertices)
				{}
				reducedGraphWithProbabilities outputGraph;
				std::vector<int> edgeCounts;
				std::size_t nUnreducedEdges;
				std::vector<int> components;
				boost::detail::depth_first_visit_restricted_impl_helper<Context::internalGraph>::stackType stack;
				std::vector<boost::default_color_type> colorMap;
				const std::vector<int>& interestVertices;
				std::vector<int> reducedInterestVertices;
			};
			void getReducedGraphNoSelfWithWeights(getReducedGraphNoSelfWithWeightsInput& input) const;
			double getRadius() const;
			subObs& operator=(subObs&&);
		protected:
			subObs(const Context& context, boost::shared_array<edgeState> state, double radius);
			subObs(subObs&& other);
			double radius;
		};
	}
}
#endif
