#include "NetworkReliabilitySubObsTree.h"
#include <graphviz/gvc.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
namespace networkReliability
{
	NetworkReliabilitySubObsTree::NetworkReliabilitySubObsTree(boost::archive::binary_iarchive& ar)
		:externalContext(NULL)
	{
		ar >> *this;
	}
	NetworkReliabilitySubObsTree::NetworkReliabilitySubObsTree(boost::archive::text_iarchive& ar)
		:externalContext(NULL)
	{
		ar >> *this;
	}
	void NetworkReliabilitySubObsTree::layout() const
	{
		int totalVertices = 0;
		for(std::vector<NetworkReliabilitySubObsCollection>::const_iterator i = levelData.begin(); i != levelData.end(); i++)
		{
			totalVertices += i->getSampleSize();
		}
		treeGraph.reset(new treeGraphType(totalVertices));
		std::size_t nLevels = levelData.size();

		//Set vertex properties, and also add edges
		treeGraphType::vertex_descriptor currentVertex = *(boost::vertices(*treeGraph).first);
		std::size_t cumulativeVertices = -levelData[0].getSampleSize();
		for(std::size_t currentLevel = 0; currentLevel < nLevels; currentLevel++)
		{
			std::size_t currentLevelSize = levelData[currentLevel].getSampleSize();
			for(std::size_t currentLevelIndex = 0; currentLevelIndex < currentLevelSize; currentLevelIndex++)
			{
				(*treeGraph)[currentVertex].potentiallyDisconnected = potentiallyDisconnected[currentLevel][currentLevelIndex];
				(*treeGraph)[currentVertex].level = currentLevel;
				(*treeGraph)[currentVertex].index = currentLevelIndex;
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
			Agraph_t* graphvizGraph = agmemread(ss.str().c_str());
				gvLayout(gvc, graphvizGraph, "dot");
					gvRender(gvc, graphvizGraph, "dot", NULL);
				gvFreeLayout(gvc, graphvizGraph);
				gvLayout(gvc, graphvizGraph, "dot");
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
	}
	const std::vector<std::vector<int > >& NetworkReliabilitySubObsTree::getPerLevelVertexIds() const
	{
		return perLevelVertexIds;
	}
	void NetworkReliabilitySubObsTree::perLevelVertexIdsFromGraph() const
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
			perLevelVertexIds[(*treeGraph)[*currentVertexIterator].level][(*treeGraph)[*currentVertexIterator].index] = *currentVertexIterator;
		}
	}
	const Context& NetworkReliabilitySubObsTree::getContext() const
	{
		if(externalContext) return *externalContext;
		return *containedContext;
	}
	NetworkReliabilitySubObsTree::NetworkReliabilitySubObsTree(Context const* externalContext, const std::vector<double>& thresholds)
		:externalContext(externalContext), thresholds(thresholds)
	{
		std::size_t nLevels = thresholds.size();
		parentData.resize(nLevels);
		potentiallyDisconnected.resize(nLevels);
		for(unsigned int i = 0; i < nLevels; i++)
		{
			NetworkReliabilitySubObsCollection currentLevelData(externalContext, thresholds[i]);
			levelData.push_back(std::move(currentLevelData));
		}
	}
	void NetworkReliabilitySubObsTree::reserve(unsigned int reservePerLevel)
	{
		std::size_t nLevels = levelData.size();
		for(std::size_t i = 0; i < nLevels; i++)
		{
			levelData[i].reserve(reservePerLevel);
			parentData[i].reserve(reservePerLevel);
			potentiallyDisconnected[i].reserve(reservePerLevel);
		}
	}
	std::size_t NetworkReliabilitySubObsTree::nLevels() const
	{
		return levelData.size();
	}
	std::size_t NetworkReliabilitySubObsTree::getSampleSize(unsigned int level) const
	{
		return levelData[level].getSampleSize();
	}
	void NetworkReliabilitySubObsTree::expand(boost::shared_array<EdgeState> state, unsigned int level, unsigned int index) const
	{
		if(level >= levelData.size())
		{
			throw std::runtime_error("Specified level does not exist in call to NetworkReliabilitySubObsTree::expand");
		}
		if(index >= levelData[level].getSampleSize())
		{
			throw std::runtime_error("Specified observation does not exist in specified level, in call to NetworkReliabilitySubObsTree::expand");
		}
		levelData[level].expand(index, state);
	}
	void NetworkReliabilitySubObsTree::add(const NetworkReliabilitySubObs& subObs, unsigned int level, int parentIndex, bool potentiallyDisconnected)
	{
		//Graph will need to be laid out again
		if(treeGraph) treeGraph.reset();
		if(level >= levelData.size())
		{
			throw std::runtime_error("Specified level does not exist in call to NetworkReliabilitySubObsTree::add");
		}
		levelData[level].add(subObs);
		parentData[level].push_back(parentIndex);
		this->potentiallyDisconnected[level].push_back(potentiallyDisconnected);
	}
	void NetworkReliabilitySubObsTree::vectorsFromGraph()
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
	const std::vector<double>& NetworkReliabilitySubObsTree::getThresholds() const
	{
		return thresholds;
	}
	const NetworkReliabilitySubObsTree::treeGraphType& NetworkReliabilitySubObsTree::getTreeGraph() const
	{
		if(!treeGraph) layout();
		return *treeGraph;
	}
}
