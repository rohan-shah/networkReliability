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
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/serialization.hpp>
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
		friend class boost::serialization::access;
		NetworkReliabilitySubObs(Context const& context, boost::archive::text_iarchive& ar);
		NetworkReliabilitySubObs(Context const& context, boost::archive::binary_iarchive& ar);
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
		void getReducedGraph(Context::internalGraph& outputGraph, int& minimumInoperative, std::vector<int>& edgeCounts, std::vector<int>& components, boost::detail::depth_first_visit_restricted_impl_helper<Context::internalGraph>::stackType& stack, std::vector<boost::default_color_type>& colorMap) const;
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
		NetworkReliabilitySubObs copyWithGeneratedObservationConditioningProb(const conditioning_type& conditioningProb) const;
		Context const& getContext() const;
		int getRadius() const;
		BOOST_SERIALIZATION_SPLIT_MEMBER()
		template<class Archive> void save(Archive& ar, const unsigned int version) const
		{
			std::string typeString = "networkReliabilitySubObs";
			ar << typeString;
			ar << boost::serialization::make_array(state.get(), context.getNEdges());
			ar << radius << minCut << couldBeDeactivated << conditioningCount << fixedInop << conditioningProb;
			typeString = "networkReliabilitySubObs_end";
			ar << typeString;
		}
		template<class Archive> void load(Archive& ar, const unsigned int version)
		{
			std::string typeString;
			ar >> typeString;
			if(typeString != "networkReliabilitySubObs")
			{
				throw std::runtime_error("Incorrect type specifier");
			}

			state.reset(new EdgeState[context.getNEdges()]);
			ar >> boost::serialization::make_array(state.get(), context.getNEdges());
			ar >> radius >> minCut >> couldBeDeactivated >> conditioningCount >> fixedInop >> conditioningProb;
			ar >> typeString;
			if(typeString != "networkReliabilitySubObs_end")
			{
				throw std::runtime_error("Incorrect type specifier");
			}
		}
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
	struct NetworkReliabilitySubObsWithContext
	{
	public:
		NetworkReliabilitySubObsWithContext(boost::archive::text_iarchive& ar)
		{
			ar >> *this;
		}
		NetworkReliabilitySubObsWithContext(boost::archive::binary_iarchive& ar)
		{
			ar >> *this;
		}
		NetworkReliabilitySubObsWithContext(NetworkReliabilitySubObs& subObs);
		friend class boost::serialization::access;
		const NetworkReliabilitySubObs& getSubObs();
	private:
		BOOST_SERIALIZATION_SPLIT_MEMBER()
		template<class Archive> void save(Archive& ar, const unsigned int version) const
		{
			std::string typeString = "networkReliabilitySubObsWithContext";
			ar << typeString;
			ar << subObs->getContext();
			ar << *subObs;
			typeString = "networkReliabilitySubObsWithContext_end";
			ar << typeString;
		}
		template<class Archive> void load(Archive& ar, const unsigned int version)
		{
			std::string typeString;
			ar >> typeString;
			if(typeString != "networkReliabilitySubObsWithContext")
			{
				throw std::runtime_error("Incorrect type specifier");
			}
			context.reset(new Context(ar));
			subObs.reset(new NetworkReliabilitySubObs(*context.get(), ar));
			ar >> typeString;
			if(typeString != "networkReliabilitySubObsWithContext_end")
			{
				throw std::runtime_error("Incorrect type specifier");
			}
		}
		boost::shared_ptr<const Context> context;
		boost::shared_ptr<NetworkReliabilitySubObs> subObs;
	};
}

#endif
