//
//=======================================================================
// Copyright 1997-2001 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//=======================================================================
//
#ifndef BOOST_GRAPH_CONNECTED_COMPONENTS_RESTRICTED_HPP
#define BOOST_GRAPH_CONNECTED_COMPONENTS_RESTRICTED_HPP

#include <boost/config.hpp>
#include "EdgeState.h"
#include "depth_first_search_restricted.hpp"
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/overloading.hpp>
#include <boost/static_assert.hpp>
#include <boost/concept/assert.hpp>
#include <boost/graph/connected_components.hpp>

namespace boost {

  template <class Graph, class ComponentMap, class ColorMap>
  inline typename property_traits<ComponentMap>::value_type
  connected_components_restricted(const Graph& g, ComponentMap c, ColorMap color, typename detail::depth_first_visit_restricted_impl_helper<Graph>::stackType& stack,
						const networkReliability::EdgeState* state
                       BOOST_GRAPH_ENABLE_IF_MODELS_PARM(Graph, vertex_list_graph_tag))
  {
    if (num_vertices(g) == 0) return 0;

    typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
    BOOST_CONCEPT_ASSERT(( WritablePropertyMapConcept<ComponentMap, Vertex> ));
    // BOOST_STATIC_ASSERT((boost::is_same<directed, undirected_tag>::value));

    typedef typename property_traits<ComponentMap>::value_type comp_type;
    // c_count initialized to "nil" (with nil represented by (max)())
    comp_type c_count((std::numeric_limits<comp_type>::max)());
    detail::components_recorder<ComponentMap> vis(c, c_count);
    depth_first_search_restricted(g, vis, color, stack, state);
    if(c_count == (std::numeric_limits<comp_type>::max)()) return 0;
    return c_count + 1;
  }

  
} // namespace boost

#ifdef BOOST_GRAPH_USE_MPI
#  include <boost/graph/distributed/connected_components.hpp>
#endif

#endif // BOOST_GRAPH_CONNECTED_COMPONENTS_HPP
