#ifndef SERIES_PARALLEL_HEADER_GUARD
#define SERIES_PARALLEL_HEADER_GUARD
namespace networkReliability
{
	template<class T> void seriesParallelReduction(T& graph, const std::vector<int>& interestVertices)
	{
		bool shouldContinue = true;
		while(shouldContinue)
		{
			shouldContinue = false;
			typename T::vertex_iterator currentVertex, end;
			boost::tie(currentVertex, end) = boost::vertices(graph);
			for(; currentVertex != end; currentVertex++)
			{
				int originalVertex = boost::get(boost::vertex_name, graph, *currentVertex);
				bool isInterestVertex = false;
				for(std::vector<int>::const_iterator i = interestVertices.begin(); i != interestVertices.end(); i++)
				{
					isInterestVertex |= (originalVertex == *i);
				}
				if(isInterestVertex) continue;
				if(boost::degree(*currentVertex, graph) == 1)
				{
					boost::clear_vertex(*currentVertex, graph);
					boost::remove_vertex(*currentVertex, graph);
					shouldContinue = true;
					break;
				}
				if(boost::degree(*currentVertex, graph) == 2)
				{
					typename T::out_edge_iterator firstDeletedEdge, secondDeletedEdge, endOutEdges;
					boost::tie(firstDeletedEdge, endOutEdges) = boost::out_edges(*currentVertex, graph);
					secondDeletedEdge = firstDeletedEdge + 1;
					typename T::vertex_descriptor firstVertex = boost::target(*firstDeletedEdge, graph), secondVertex = boost::target(*secondDeletedEdge, graph);
					if(firstVertex == secondVertex) 
					{
						boost::clear_vertex(*currentVertex, graph);
						boost::remove_vertex(*currentVertex, graph);
						shouldContinue = true;
						break;
					}
				
					mpfr_class firstOp = boost::get(boost::edge_op_probability, graph, *firstDeletedEdge);
					mpfr_class firstInop = boost::get(boost::edge_inop_probability, graph, *firstDeletedEdge);
					mpfr_class secondOp = boost::get(boost::edge_op_probability, graph, *secondDeletedEdge);
					mpfr_class secondInop = boost::get(boost::edge_inop_probability, graph, *secondDeletedEdge);
					
					mpfr_class newOp = firstOp*secondOp;
					mpfr_class newInop = 1 - newOp;
					
					bool edgeAlreadyExists = false;
					typename T::out_edge_iterator currentExistingEdge;
					boost::tie(currentExistingEdge, endOutEdges) = boost::out_edges(firstVertex, graph);
					for(;currentExistingEdge != endOutEdges; currentExistingEdge++)
					{
						if(boost::target(*currentExistingEdge, graph) == secondVertex)
						{
							newInop *= boost::get(boost::edge_inop_probability, graph, *currentExistingEdge);
							newOp = 1 - newInop;
							boost::put(boost::edge_op_probability, graph, *currentExistingEdge, newOp);
							boost::put(boost::edge_inop_probability, graph, *currentExistingEdge, newInop);
			
							edgeAlreadyExists = true;
							break;
						}
					}
					if(!edgeAlreadyExists)
					{
						typename T::edge_descriptor newEdge = boost::add_edge(firstVertex, secondVertex, graph).first;
						boost::put(boost::edge_op_probability, graph, newEdge, newOp);
						boost::put(boost::edge_inop_probability, graph, newEdge, newInop);
					}

					boost::clear_vertex(*currentVertex, graph);
					boost::remove_vertex(*currentVertex, graph);
					shouldContinue = true;
					break;
				}
			}
		}
	}
}
#endif
