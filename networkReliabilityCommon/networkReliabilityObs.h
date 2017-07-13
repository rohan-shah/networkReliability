#ifndef NETWORK_RELIABILITY_OBS_HEADER_GUARD
#define NETWORK_RELIABILITY_OBS_HEADER_GUARD
#include "context.h"
#include "edgeState.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/shared_array.hpp>
#include <boost/noncopyable.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
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
	class NetworkReliabilityObs : public boost::noncopyable
	{
	public:
		friend class boost::serialization::access;
		NetworkReliabilityObs(context const& contextObj, boost::mt19937& randomSource);
		NetworkReliabilityObs(context const& contextObj, boost::shared_array<edgeState> state);
		NetworkReliabilityObs(context const& contextObj, boost::archive::binary_iarchive& archive);
		NetworkReliabilityObs(context const& contextObj, boost::archive::text_iarchive& archive);
		NetworkReliabilityObs(NetworkReliabilityObs&& other);
		const edgeState* getState() const;
		NetworkReliabilityObs& operator=(const NetworkReliabilityObs& other);
		const context& getContext() const;
		//The reduced graph code
		void getReducedGraph(context::internalGraph& outputGraph, int& minimumInoperative, std::vector<int>& edgeCounts, std::vector<int>& components, boost::detail::depth_first_visit_restricted_impl_helper<context::internalGraph>::stackType& stack, std::vector<boost::default_color_type>& colorMap) const;
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
			boost::detail::depth_first_visit_restricted_impl_helper<context::internalGraph>::stackType stack;
			std::vector<boost::default_color_type> colorMap;
			const std::vector<int>& interestVertices;
			std::vector<int> reducedInterestVertices;
		};
		void getReducedGraphNoSelfWithWeights(getReducedGraphNoSelfWithWeightsInput& input) const;

	protected:
		context const& contextObj;
		boost::shared_array<edgeState> state;
	private:
		BOOST_SERIALIZATION_SPLIT_MEMBER()
		template<class Archive> void save(Archive& ar, const unsigned int version) const
		{
			std::string typeString = "networkReliabilityObs";
			ar << typeString;
			ar << boost::serialization::make_array(state.get(), contextObj.getNEdges());
			typeString = "networkReliabilityObs_end";
			ar << typeString;
		}
		template<class Archive> void load(Archive& ar, const unsigned int version)
		{
			std::string typeString;
			ar >> typeString;
			if(typeString != "networkReliabilityObs")
			{
				throw std::runtime_error("Incorrect type specifier");
			}
			state.reset(new edgeState[contextObj.getNEdges()]);
			ar >> boost::serialization::make_array(state.get(), contextObj.getNEdges());
			ar >> typeString;
			if(typeString != "networkReliabilityObs_end")
			{
				throw std::runtime_error("Incorrect type specifier");
			}
		}
	};
	struct NetworkReliabilityObsWithContext
	{
	public:
		NetworkReliabilityObsWithContext(boost::archive::text_iarchive& ar)
		{
			ar >> *this;
		}
		NetworkReliabilityObsWithContext(boost::archive::binary_iarchive& ar)
		{
			ar >> *this;
		}
		NetworkReliabilityObsWithContext(NetworkReliabilityObs& obs);
		friend class boost::serialization::access;
		const NetworkReliabilityObs& getObs() const;
		const context& getContext() const;
	private:
		BOOST_SERIALIZATION_SPLIT_MEMBER()
		template<class Archive> void save(Archive& ar, const unsigned int version) const
		{
			std::string typeString = "networkReliabilityObsWithContext";
			ar << typeString;
			ar << obs->getContext();
			ar << *obs;
			typeString = "networkReliabilityObsWithContext_end";
			ar << typeString;
		}
		template<class Archive> void load(Archive& ar, const unsigned int version)
		{
			std::string typeString;
			ar >> typeString;
			if(typeString != "networkReliabilityObsWithContext")
			{
				throw std::runtime_error("Incorrect type specifier");
			}
			contextPtr.reset(new context(ar));
			obs.reset(new NetworkReliabilityObs(*contextPtr.get(), ar));
			ar >> typeString;
			if(typeString != "networkReliabilityObsWithContext_end")
			{
				throw std::runtime_error("Incorrect type specifier");
			}
		}
		boost::shared_ptr<const context> contextPtr;
		boost::shared_ptr<NetworkReliabilityObs> obs;
	};

}
#endif
