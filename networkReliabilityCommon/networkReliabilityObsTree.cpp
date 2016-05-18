#include "networkReliabilityObsTree.h"
#ifdef HAS_GRAPHVIZ
	#ifdef _MSC_VER
		#define WIN32_DLL
	#endif
	#include <graphviz/gvc.h>
	#undef WIN32_DLL
#endif
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
namespace networkReliability
{
	NetworkReliabilityObsTree::NetworkReliabilityObsTree(boost::archive::binary_iarchive& ar)
		:externalContext(NULL)
	{
		ar >> *this;
	}
	NetworkReliabilityObsTree::NetworkReliabilityObsTree(boost::archive::text_iarchive& ar)
		:externalContext(NULL)
	{
		ar >> *this;
	}
	bool NetworkReliabilityObsTree::layout() const
	{
#ifdef HAS_GRAPHVIZ
		std::size_t totalVertices = 0;
		for(std::vector<NetworkReliabilityObsCollection>::const_iterator i = levelData.begin(); i != levelData.end(); i++)
		{
			totalVertices += i->getSampleSize();
		}
		treeGraph.reset(new treeGraphType(totalVertices));
		std::size_t nLevels = levelData.size();

		//Set vertex properties, and also add edges
		treeGraphType::vertex_descriptor currentVertex = *(boost::vertices(*treeGraph).first);
		std::ptrdiff_t cumulativeVertices = -(std::ptrdiff_t)levelData[0].getSampleSize();
		for(std::size_t currentLevel = 0; currentLevel < nLevels; currentLevel++)
		{
			std::size_t currentLevelSize = levelData[currentLevel].getSampleSize();
			for(std::size_t currentLevelIndex = 0; currentLevelIndex < currentLevelSize; currentLevelIndex++)
			{
				(*treeGraph)[currentVertex].potentiallyDisconnected = potentiallyDisconnected[currentLevel][currentLevelIndex];
				(*treeGraph)[currentVertex].level = (int)currentLevel;
				(*treeGraph)[currentVertex].index = (int)currentLevelIndex;
				if(parentData[currentLevel][currentLevelIndex] >= 0)
				{
					boost::add_edge(cumulativeVertices + parentData[currentLevel][currentLevelIndex], currentVertex, *treeGraph);
				}
				currentVertex++;
			}
			cumulativeVertices += currentLevelSize;
		}


		//Construct a graphviz string representation
		std::stringstream ss;
		class propertyWriter
		{
		public:
			propertyWriter(treeGraphType& t)
				:t(t)
			{}
			void operator()(std::ostream& out, const treeGraphType::vertex_descriptor& v) const
			{
				out << "[potDiscon=\"" << t[v].potentiallyDisconnected <<"\" level=\"" << t[v].level << "\" index=\"" << t[v].index << "\"]";
			}
		private:
			treeGraphType& t;
		};
		propertyWriter p(*treeGraph);
		boost::write_graphviz(ss, *treeGraph, p);
		
		//Now have graphviz create a graph from that string
		GVC_t* gvc = gvContext();
			//Older versions of graphviz appear to use just a char*
			Agraph_t* graphvizGraph = agmemread((char*)ss.str().c_str());
				int graphVizResult = gvLayout(gvc, graphvizGraph, "dot");
				if(graphVizResult == -1) return false;
					graphVizResult = gvRender(gvc, graphvizGraph, "dot", NULL);
					if(graphVizResult == -1) return false;
				gvFreeLayout(gvc, graphvizGraph);
				graphVizResult = gvLayout(gvc, graphvizGraph, "dot");
				if(graphVizResult == -1) return false;
					char* laidOutDot;
					unsigned int laidOutDotLength;
					gvRenderData(gvc, graphvizGraph, "dot", &laidOutDot, &laidOutDotLength);
						//Set up property maps for reading back in laid out graph
						boost::dynamic_properties dynamicProperties(boost::ignore_other_properties);
						boost::vector_property_map<std::string> posProperty;
						boost::property_map<treeGraphType, int vertexProperty::*>::type levelProperty(boost::get(&vertexProperty::level, *treeGraph));
						boost::property_map<treeGraphType, int vertexProperty::*>::type indexProperty(boost::get(&vertexProperty::index, *treeGraph));
						boost::property_map<treeGraphType, bool vertexProperty::*>::type potentiallyDisconnectedProperty(boost::get(&vertexProperty::potentiallyDisconnected, *treeGraph));
						dynamicProperties.property("pos", posProperty);
						dynamicProperties.property("level", levelProperty);
						dynamicProperties.property("index", indexProperty);
						dynamicProperties.property("potDiscon", potentiallyDisconnectedProperty);
						//actually read graph back into boos
						treeGraph->clear();
						boost::read_graphviz(laidOutDot, laidOutDot+laidOutDotLength, *treeGraph, dynamicProperties);
						//Now split pos up into x and y from pos, and put those into the graph
						{
							treeGraphType::vertex_iterator current, end;
							boost::tie(current, end) = boost::vertices(*treeGraph);
							std::vector<std::string> parts;
							for(;current != end; current++)
							{
								parts.clear();
								boost::algorithm::split(parts, posProperty[*current], boost::algorithm::is_any_of(","), boost::token_compress_on);
								double x = boost::lexical_cast<double>(parts[0]);
								double y = boost::lexical_cast<double>(parts[1]);
								(*treeGraph)[*current].x = x;
								(*treeGraph)[*current].y = y;
							}
						}
					gvFreeRenderData(laidOutDot);
				gvFreeLayout(gvc, graphvizGraph);
			agclose(graphvizGraph);
		gvFreeContext(gvc);
		perLevelVertexIdsFromGraph();
		return true;
#else 
		return true;
#endif
	}
	const std::vector<std::vector<int > >& NetworkReliabilityObsTree::getPerLevelVertexIds() const
	{
		return perLevelVertexIds;
	}
	void NetworkReliabilityObsTree::perLevelVertexIdsFromGraph() const
	{
		std::size_t nLevels = levelData.size();
		perLevelVertexIds.clear();
		perLevelVertexIds.resize(nLevels);
		for(std::size_t level = 0; level < nLevels; level++)
		{
			perLevelVertexIds[level].resize(levelData[level].getSampleSize());
		}
		treeGraphType::vertex_iterator currentVertexIterator, endVertexIterator;
		boost::tie(currentVertexIterator, endVertexIterator) = boost::vertices(*treeGraph);
		for(; currentVertexIterator != endVertexIterator; currentVertexIterator++)
		{
			perLevelVertexIds[(*treeGraph)[*currentVertexIterator].level][(*treeGraph)[*currentVertexIterator].index] = (int)*currentVertexIterator;
		}
	}
	const context& NetworkReliabilityObsTree::getContext() const
	{
		if(externalContext) return *externalContext;
		return *containedContext;
	}
	NetworkReliabilityObsTree::NetworkReliabilityObsTree(context const* externalContext, const std::vector<double>& thresholds)
		:externalContext(externalContext), thresholds(thresholds)
	{
		std::size_t nLevels = thresholds.size();
		parentData.resize(nLevels);
		potentiallyDisconnected.resize(nLevels);
		for(unsigned int i = 0; i < nLevels; i++)
		{
			NetworkReliabilityObsCollection currentLevelData(externalContext, thresholds[i]);
			levelData.push_back(std::move(currentLevelData));
		}
	}
	void NetworkReliabilityObsTree::reserve(std::size_t reservePerLevel)
	{
		std::size_t nLevels = levelData.size();
		for(std::size_t i = 0; i < nLevels; i++)
		{
			levelData[i].reserve(reservePerLevel);
			parentData[i].reserve(reservePerLevel);
			potentiallyDisconnected[i].reserve(reservePerLevel);
		}
	}
	std::size_t NetworkReliabilityObsTree::nLevels() const
	{
		return levelData.size();
	}
	std::size_t NetworkReliabilityObsTree::getSampleSize(unsigned int level) const
	{
		return levelData[level].getSampleSize();
	}
	void NetworkReliabilityObsTree::expand(boost::shared_array<edgeState> state, unsigned int level, unsigned int index) const
	{
		if(level >= levelData.size())
		{
			throw std::runtime_error("Specified level does not exist in call to NetworkReliabilityObsTree::expand");
		}
		if(index >= levelData[level].getSampleSize())
		{
			throw std::runtime_error("Specified observation does not exist in specified level, in call to NetworkReliabilityObsTree::expand");
		}
		levelData[level].expand(index, state);
	}
	void NetworkReliabilityObsTree::add(const NetworkReliabilityObs& obs, unsigned int level, int parentIndex, bool potentiallyDisconnected)
	{
		//Graph will need to be laid out again
		if(treeGraph) treeGraph.reset();
		if(level >= levelData.size())
		{
			throw std::runtime_error("Specified level does not exist in call to NetworkReliabilityObsTree::add");
		}
		levelData[level].add(obs);
		parentData[level].push_back(parentIndex);
		this->potentiallyDisconnected[level].push_back(potentiallyDisconnected);
	}
	void NetworkReliabilityObsTree::vectorsFromGraph()
	{
		assert(treeGraph);
		treeGraphType::vertex_iterator current, end;
		boost::tie(current, end) = boost::vertices(*treeGraph);

		std::size_t nLevels = levelData.size();
		potentiallyDisconnected.clear();
		parentData.clear();
		potentiallyDisconnected.resize(nLevels);
		parentData.resize(nLevels);

		for(std::size_t level = 0; level < nLevels; level++)
		{
			parentData[level].resize(levelData[level].getSampleSize());
			potentiallyDisconnected[level].resize(levelData[level].getSampleSize());
		}
		for(; current != end; current++)
		{
			//This one we have to actually get the index of the parent vertex
			treeGraphType::in_edge_iterator current_in, end_in;
			boost::tie(current_in, end_in) = boost::in_edges(*current, *treeGraph);
			int level = (*treeGraph)[*current].level;
			int index = (*treeGraph)[*current].index;
			if(current_in == end_in) parentData[level][index] = -1;
			else parentData[level][index] = (*treeGraph)[current_in->m_source].index;
			//This one is just a lookup
			potentiallyDisconnected[level][(*treeGraph)[*current].index] = (*treeGraph)[*current].potentiallyDisconnected;
		}
		perLevelVertexIdsFromGraph();
	}
	const std::vector<double>& NetworkReliabilityObsTree::getThresholds() const
	{
		return thresholds;
	}
	const NetworkReliabilityObsTree::treeGraphType& NetworkReliabilityObsTree::getTreeGraph() const
	{
		if(!treeGraph) 
		{
			//This returns true in the case of graphviz not available
			bool result = layout();
			if(!result)
			{
				//This path is taken if graphviz is available but the layout fails
				throw std::runtime_error("Call to NetworkReliabilityObsTree::layout failed");
			}
		}
		if(!treeGraph)
		{
			throw std::runtime_error("Call to NetworkReliabilityObsTree::getTreeGraph failed, possibly because graphviz was not available");
		}
		return *treeGraph;
	}
}
