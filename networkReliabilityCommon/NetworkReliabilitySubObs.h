#ifndef NETWORK_RELIABILITY_SUB_OBS_HEADER_GUARD
#define NETWORK_RELIABILITY_SUB_OBS_HEADER_GUARD
#include <boost/random/mersenne_twister.hpp>
#include "Context.h"
#include "serializeGMP.hpp"
#include "EdgeState.h"
#include <boost/shared_array.hpp>
#include "includeMPFR.h"
#include "graphAlgorithms.h"
#include <boost/serialization/access.hpp>
#include <boost/serialization/string.hpp>
namespace boost
{
#define BOOST_INSTALL_PROPERTY(KIND, NAME) \
  template <> struct property_kind<KIND##_##NAME##_t> { \
    typedef KIND##_property_tag type; \
  }

#define BOOST_DEF_PROPERTY(KIND, NAME) \
  enum KIND##_##NAME##_t { KIND##_##NAME }; \
  BOOST_INSTALL_PROPERTY(KIND, NAME)
	BOOST_DEF_PROPERTY(edge, inop_probability);
	BOOST_DEF_PROPERTY(edge, op_probability);
}
namespace networkReliability
{
	class NetworkReliabilityObs;
	class NetworkReliabilitySubObs : public boost::noncopyable
	{
	public:
		template<class Archive> friend void readNetworkReliabilitySubObs(Archive& ar, boost::shared_ptr<const Context>& context, boost::shared_ptr<NetworkReliabilitySubObs>& subObs);
		template<class Archive> friend void writeNetworkReliabilitySubObs(Archive& ar, NetworkReliabilitySubObs& subObs);
		typedef mpfr_class conditioning_type;
		NetworkReliabilitySubObs(NetworkReliabilitySubObs&& other);
		NetworkReliabilitySubObs(Context const& context, boost::shared_array<EdgeState> state, int radius, int conditioningCount, conditioning_type conditiniongProb);
		NetworkReliabilityObs getObservation(boost::mt19937& randomSource) const;
		const EdgeState* getState() const;
		int getMinCut() const;
		NetworkReliabilitySubObs& operator=(NetworkReliabilitySubObs&& other);
		int getConditioningCount() const;
		const conditioning_type& getConditioningProb() const;
		const conditioning_type& getGeneratedObservationConditioningProb() const;
		void getRadius1ReducedGraph(Context::internalGraph& outputGraph, int& minimumInoperative, std::vector<int>& edgeCounts, std::vector<int>& components, boost::detail::depth_first_visit_restricted_impl_helper<Context::internalGraph>::stackType& stack, std::vector<boost::default_color_type>& colorMap) const;
		typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::property<boost::vertex_name_t, int>, boost::property<boost::edge_index_t, int, boost::property<boost::edge_inop_probability_t, mpfr_class, boost::property<boost::edge_op_probability_t, mpfr_class> > > > reducedGraphWithProbabilities;
		struct getRadius1ReducedGraphNoSelfWithWeightsInput
		{
		public:
			getRadius1ReducedGraphNoSelfWithWeightsInput(const std::vector<int>& interestVertices)
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
		void getRadius1ReducedGraphNoSelfWithWeights(getRadius1ReducedGraphNoSelfWithWeightsInput& input) const;
		NetworkReliabilitySubObs copyWithGeneratedObservationConditioningProb(const conditioning_type& conditioningProb) const;
		Context const& getContext() const;
		int getRadius() const;
	private:
		NetworkReliabilitySubObs(Context const& context);
		Context const& context;
		boost::shared_array<EdgeState> state;
		int radius;
		int minCut;
		std::vector<int> couldBeDeactivated;
		int conditioningCount;
		int fixedInop;
		conditioning_type conditioningProb, generatedObservationConditioningProb;
		int generatedObservationConditioningCount;
	};
	template<class Archive> void readNetworkReliabilitySubObs(Archive& ar, boost::shared_ptr<const Context>& context, boost::shared_ptr<NetworkReliabilitySubObs>& subObs)
	{
		std::string typeString;
		ar >> typeString;
		if(typeString != "networkReliabilitySubObs")
		{
			throw std::runtime_error("File did not start with correct type specifier");
		}
		boost::shared_ptr<Context> nonConstContext(new Context());
		ar >> *nonConstContext;
		context = nonConstContext;

		subObs.reset(new NetworkReliabilitySubObs(*context));
		NetworkReliabilitySubObs& subObsRef = *subObs;
		subObsRef.state.reset(new EdgeState[context->getNEdges()]);
		ar >> boost::serialization::make_array(subObsRef.state.get(), context->getNEdges());
		ar >> subObsRef.radius >> subObsRef.minCut >> subObsRef.couldBeDeactivated >> subObsRef.conditioningCount >> subObsRef.fixedInop >> subObsRef.conditioningProb;
	}
	template<class Archive> void writeNetworkReliabilitySubObs(Archive& ar, NetworkReliabilitySubObs& subObs)
	{
		std::string typeString = "networkReliabilitySubObs";
		ar << typeString;
		ar << subObs.getContext();
		ar << boost::serialization::make_array(subObs.state.get(), subObs.getContext().getNEdges());
		ar << subObs.radius << subObs.minCut << subObs.couldBeDeactivated << subObs.conditioningCount << subObs.fixedInop << subObs.conditioningProb;
	}
}
#endif
