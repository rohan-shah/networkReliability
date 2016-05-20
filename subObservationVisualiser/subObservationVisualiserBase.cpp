//I want assert always enabled in this particular file
#ifdef NDEBUG
#undef NDEBUG
#include <assert.h>
#define NDEBUG
#endif

#include "subObservationVisualiserBase.h"
#include <QEvent>
#include <QKeyEvent>
#include <QTimer>
#include "ZoomGraphicsView.h"
#include <QGraphicsRectItem>
#include "graphAlgorithms.h"
#include <boost/lexical_cast.hpp>
#include <QGraphicsSceneMouseEvent>
namespace networkReliability
{
	bool sortByFirst(context::vertexPosition const& first, context::vertexPosition const& second)
	{
		return first.first < second.first;
	}
	bool sortBySecond(context::vertexPosition const& first, context::vertexPosition const& second)
	{
		return first.second < second.second;
	}
	void subObservationVisualiserBase::setReduced(bool reduced)
	{
		if(reduced != this->reduced)
		{
			this->reduced = reduced;
			if(!this->reduced) highlightedReducedComponent = -1;
			updateGraphics();
		}
	}
	void subObservationVisualiserBase::switchReduced()
	{
		reduced = !reduced;
		if(!reduced) highlightedReducedComponent = -1;
		updateGraphics();
	}
	void subObservationVisualiserBase::updateReducedGraphData(const NetworkReliabilityObs& subObs)
	{
		boost::detail::depth_first_visit_restricted_impl_helper<context::internalGraph>::stackType stack;
		std::vector<boost::default_color_type> colorMap;
		nUnreducedComponents = countComponents(contextObj, subObs.getState(), reducedComponents, stack, colorMap);
		subObs.getReducedGraphNoSelfWithWeights(reducedGraphData);
	}
	void subObservationVisualiserBase::initialiseControls()
	{
		graphicsScene = new QGraphicsScene();
		graphicsScene->installEventFilter(this);
		graphicsScene->setItemIndexMethod(QGraphicsScene::NoIndex);

		graphicsView = new ZoomGraphicsView(graphicsScene);
		graphicsView->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
		graphicsView->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
		graphicsView->viewport()->installEventFilter(this);

		const std::vector<context::vertexPosition>& vertexPositions = contextObj.getVertexPositions();
		minX = std::min_element(vertexPositions.begin(), vertexPositions.end(), sortByFirst)->first - pointSize;
		maxX = std::max_element(vertexPositions.begin(), vertexPositions.end(), sortByFirst)->first + pointSize;
		minY = std::min_element(vertexPositions.begin(), vertexPositions.end(), sortBySecond)->second - pointSize;
		maxY = std::max_element(vertexPositions.begin(), vertexPositions.end(), sortBySecond)->second + pointSize;

		layout = new QHBoxLayout;
		layout->addWidget(graphicsView, 1);
		layout->setContentsMargins(0,0,0,0);
		setLayout(layout);
	}
	subObservationVisualiserBase::subObservationVisualiserBase(const context& contextObj, float pointSize)
		:pointSize(pointSize), highlightedReducedComponent(-1), reducedGraphData(contextObj.getInterestVertices()), reduced(false), contextObj(contextObj), reducedPointsItem(NULL), reducedLinesItem(NULL), unreducedPointsItem(NULL), unreducedLinesItem(NULL), highlightedGroupItem(NULL), backgroundItem(NULL)
	{
		initialiseControls();

		QPen pen(Qt::NoPen);
		QColor grey("grey");
		grey.setAlphaF(0.5);
		QBrush brush;
		brush.setColor(grey);
		brush.setStyle(Qt::SolidPattern);
		backgroundItem = new QGraphicsRectItem(minX, minY, maxX - minX, maxY - minY, NULL);
		backgroundItem->setPen(pen);
		backgroundItem->setBrush(brush);
		backgroundItem->setZValue(-1);
		graphicsScene->addItem(backgroundItem);

		unreducedPointsItem = new QGraphicsItemGroup(NULL);
		constructUnreducedPoints();
	}
	void subObservationVisualiserBase::setObservation(const NetworkReliabilityObs& subObs)
	{
		highlightedReducedComponent = -1;
		updateReducedGraphData(subObs);
		constructGraphics(subObs);
		updateGraphics();
	}
	subObservationVisualiserBase::~subObservationVisualiserBase()
	{
	}
	void subObservationVisualiserBase::updateGraphics()
	{
		QList<QGraphicsItem*> allItems = graphicsScene->items();
		for(QList<QGraphicsItem*>::iterator i = allItems.begin(); i != allItems.end(); i++) 
		{
			if(!(*i)->parentItem()) graphicsScene->removeItem(*i);
		}
		graphicsScene->addItem(backgroundItem);
		if(reduced)
		{
			graphicsScene->addItem(reducedPointsItem);
			graphicsScene->addItem(reducedLinesItem);
		}
		else
		{
			graphicsScene->addItem(unreducedPointsItem);
			graphicsScene->addItem(unreducedLinesItem);
		}
	}
	void subObservationVisualiserBase::constructGraphics(const NetworkReliabilityObs& subObs)
	{
		QList<QGraphicsItem*> allItems = graphicsScene->items();
		for(QList<QGraphicsItem*>::iterator i = allItems.begin(); i != allItems.end(); i++) 
		{
			if(!(*i)->parentItem()) graphicsScene->removeItem(*i);
		}

		delete reducedPointsItem;
		delete reducedLinesItem;
		delete unreducedLinesItem;
		delete highlightedGroupItem;

		reducedPointsItem = new QGraphicsItemGroup(NULL);
		reducedLinesItem = new QGraphicsItemGroup(NULL);
		unreducedLinesItem = new QGraphicsItemGroup(NULL);
		highlightedGroupItem = new QGraphicsItemGroup(NULL);

		constructReducedLines(subObs);
		constructUnreducedLines(subObs);
		constructReducedPoints();
	}
	void subObservationVisualiserBase::constructUnreducedPoints()
	{
		assert(unreducedPointsItem);
		std::size_t nVertices = boost::num_vertices(contextObj.getGraph());

		const std::vector<context::vertexPosition>& vertexPositions = contextObj.getVertexPositions();
		const std::vector<int>& interestVertices = contextObj.getInterestVertices();

		QPen blackPen(QColor("black"));
		blackPen.setStyle(Qt::NoPen);
		QBrush blackBrush(QColor("black"));

		QPen redPen(QColor("red"));
		redPen.setStyle(Qt::NoPen);
		QBrush redBrush(QColor("red"));
		
		QPen greenPen(QColor("green"));
		greenPen.setStyle(Qt::NoPen);
		QBrush greenBrush(QColor("green"));

		for(std::size_t vertexCounter = 0; vertexCounter < nVertices; vertexCounter++)
		{
			context::vertexPosition currentPosition = vertexPositions[vertexCounter];
			float x = currentPosition.first;
			float y = currentPosition.second;
			if(highlightedReducedComponent != -1 && reducedComponents[vertexCounter] == highlightedReducedComponent)
			{
				QGraphicsEllipseItem* newItem = new QGraphicsEllipseItem(x - pointSize, y - pointSize, 2*pointSize, 2*pointSize, unreducedPointsItem);
				newItem->setBrush(greenBrush);
				newItem->setPen(greenPen);
			}
			if(interestVertices.end() == std::find(interestVertices.begin(), interestVertices.end(), vertexCounter))
			{
				QGraphicsEllipseItem* newItem = new QGraphicsEllipseItem(x - pointSize/2, y - pointSize/2, pointSize, pointSize, unreducedPointsItem);
				newItem->setBrush(blackBrush);
				newItem->setPen(blackPen);
			}
			else
			{
				QGraphicsEllipseItem* newItem = new QGraphicsEllipseItem(x - pointSize/2, y - pointSize/2, pointSize, pointSize, unreducedPointsItem);
				newItem->setBrush(redBrush);
				newItem->setPen(redPen);
			}
		}
	}
	void subObservationVisualiserBase::constructUnreducedLines(const NetworkReliabilityObs& subObs)
	{
		assert(unreducedLinesItem);
		const EdgeState* state = subObs.getState();
		const context::internalGraph& graph = contextObj.getGraph();
		const std::vector<context::vertexPosition>& vertexPositions = contextObj.getVertexPositions();

		QPen pen(QColor("black"));
		pen.setWidthF(pointSize/10);

		QPen redPen(QColor("red"));
		redPen.setWidthF(pointSize/10);

		QVector<qreal> dashPattern;
		dashPattern.push_back(8);
		dashPattern.push_back(8);

		QPen dashedPen(QColor("black"));
		dashedPen.setDashPattern(dashPattern);

		QPen dashedRedPen(QColor("red"));
		dashedRedPen.setDashPattern(dashPattern);

		context::internalGraph::edge_iterator start, end;
		boost::tie(start, end) = boost::edges(graph);

		boost::property_map<context::internalGraph, boost::edge_index_t>::const_type edgeIndexMap = boost::get(boost::edge_index, graph);

		while(start != end)
		{
			context::vertexPosition sourcePosition = vertexPositions[start->m_source], targetPosition = vertexPositions[start->m_target];
			if(state[boost::get(boost::edge_index, graph, *start)] & FIXED_OP)
			{
				QGraphicsLineItem* newItem = new QGraphicsLineItem(sourcePosition.first, sourcePosition.second, targetPosition.first, targetPosition.second, unreducedLinesItem);
				newItem->setPen(redPen);
			}
			else if(state[boost::get(boost::edge_index, graph, *start)] & UNFIXED_OP)
			{
				QGraphicsLineItem* newItem = new QGraphicsLineItem(sourcePosition.first, sourcePosition.second, targetPosition.first, targetPosition.second, unreducedLinesItem);
				newItem->setPen(pen);
			}
			else if(state[boost::get(boost::edge_index, graph, *start)] & UNFIXED_INOP)
			{
				QGraphicsLineItem* newItem = new QGraphicsLineItem(sourcePosition.first, sourcePosition.second, targetPosition.first, targetPosition.second, unreducedLinesItem);
				newItem->setPen(dashedPen);
			}
			else
			{
				QGraphicsLineItem* newItem = new QGraphicsLineItem(sourcePosition.first, sourcePosition.second, targetPosition.first, targetPosition.second, unreducedLinesItem);
				newItem->setPen(dashedRedPen);
			}
			QGraphicsSimpleTextItem* textItem = new QGraphicsSimpleTextItem(QString::fromStdString(boost::lexical_cast<std::string>(boost::get(edgeIndexMap, *start))), unreducedLinesItem); 
			textItem->setPos((sourcePosition.first + targetPosition.first)/2, (sourcePosition.second + targetPosition.second)/2);
			start++;
		}
	}
	bool subObservationVisualiserBase::eventFilter(QObject* object, QEvent *event)
	{
		if(event->type() == QEvent::GraphicsSceneMouseMove && object == graphicsScene)
		{
			QGraphicsSceneMouseEvent* mouseEvent = static_cast<QGraphicsSceneMouseEvent*>(event);
			QPointF position = mouseEvent->scenePos();
			emit positionChanged(position.x(), position.y());

			const context::internalGraph& graph = contextObj.getGraph();
			std::size_t nVertices = boost::num_vertices(graph);
			const std::vector<context::vertexPosition>& vertexPositions = contextObj.getVertexPositions();

			if(!reduced)
			{
				std::vector<double> distances(vertexPositions.size());
				for(std::size_t i = 0; i < nVertices; i++)
				{
					distances[i] = (position.x() - vertexPositions[i].first) * (position.x() - vertexPositions[i].first) + (position.y() - vertexPositions[i].second) * (position.y() - vertexPositions[i].second);
				}
				int closestVertex = std::distance(distances.begin(), std::min_element(distances.begin(), distances.end()));
				if(highlightedReducedComponent != reducedComponents[closestVertex])
				{
					highlightedReducedComponent = reducedComponents[closestVertex];
					updateGraphics();
				}
			}
			return true;
		}
		else if(event->type() == QEvent::KeyPress)
		{
			QKeyEvent* keyEvent = static_cast<QKeyEvent*>(event);
			if(keyEvent->key() == Qt::Key_Left && keyEvent->modifiers() & Qt::ShiftModifier)
			{
				emit observationLeft();
				return true;
			}
			else if(keyEvent->key() == Qt::Key_Right && keyEvent->modifiers() & Qt::ShiftModifier)
			{
				emit observationRight();
				return true;
			}
			else if(keyEvent->key() == Qt::Key_Up && keyEvent->modifiers() & Qt::ShiftModifier)
			{
				emit observationDown();
				return true;
			}
			else if(keyEvent->key() == Qt::Key_Down && keyEvent->modifiers() & Qt::ShiftModifier)
			{
				emit observationUp();
				return true;
			}
		}
		else if(event->type() == QEvent::Leave && object == graphicsScene)
		{
			highlightedReducedComponent = -1;
			updateGraphics();
			return true;
		}
		return false;
	}
	void subObservationVisualiserBase::constructReducedPoints()
	{
		assert(reducedPointsItem);
		const context::internalGraph& unreducedGraph = contextObj.getGraph();
		const std::vector<context::vertexPosition>& vertexPositions = contextObj.getVertexPositions();

		std::vector<bool> stillPresentInReduced(nUnreducedComponents, false);
		NetworkReliabilityObs::reducedGraphWithProbabilities::vertex_iterator currentReducedVertex, endReducedVertex;
		boost::tie(currentReducedVertex, endReducedVertex) = boost::vertices(reducedGraphData.outputGraph);
		for(;currentReducedVertex != endReducedVertex; currentReducedVertex++)
		{
			stillPresentInReduced[boost::get(boost::vertex_name, reducedGraphData.outputGraph, *currentReducedVertex)] = true;
		}
		
		QPen blackPen(QColor("black"));
		blackPen.setStyle(Qt::NoPen);
		QBrush blackBrush(QColor("black"));

		context::internalGraph::vertex_iterator currentUnreducedVertex, endUnreducedVertex;
		boost::tie(currentUnreducedVertex, endUnreducedVertex) = boost::vertices(unreducedGraph);
		for(;currentUnreducedVertex != endUnreducedVertex; currentUnreducedVertex++)
		{
			if(stillPresentInReduced[reducedComponents[*currentUnreducedVertex]])
			{
				context::vertexPosition currentPosition = vertexPositions[*currentUnreducedVertex];
				float x = currentPosition.first;
				float y = currentPosition.second;

				QGraphicsEllipseItem* newItem = new QGraphicsEllipseItem(x - pointSize/2, y - pointSize/2, pointSize, pointSize, reducedPointsItem);
				newItem->setBrush(blackBrush);
				newItem->setPen(blackPen);
			}
		}
	}
	void subObservationVisualiserBase::constructReducedLines(const NetworkReliabilityObs& subObs)
	{
		assert(reducedLinesItem);
		const context::internalGraph& unreducedGraph = contextObj.getGraph();
		const std::vector<Context::vertexPosition>& vertexPositions = contextObj.getVertexPositions();

		std::vector<bool> stillPresentInReduced(nUnreducedComponents, false);
		NetworkReliabilityObs::reducedGraphWithProbabilities::vertex_iterator currentReducedVertex, endReducedVertex;
		boost::tie(currentReducedVertex, endReducedVertex) = boost::vertices(reducedGraphData.outputGraph);
		for(;currentReducedVertex != endReducedVertex; currentReducedVertex++)
		{
			stillPresentInReduced[boost::get(boost::vertex_name, reducedGraphData.outputGraph, *currentReducedVertex)] = true;
		}
		
		QPen blackPen(QColor("black"));

		Context::internalGraph::edge_iterator currentEdge, lastEdge;
		boost::tie(currentEdge, lastEdge) = boost::edges(unreducedGraph);
		const EdgeState* state = subObs.getState();
		for(; currentEdge != lastEdge; currentEdge++)
		{
			int sourceVertex = (int)boost::source(*currentEdge, unreducedGraph);
			int targetVertex = (int)boost::target(*currentEdge, unreducedGraph);
			if(stillPresentInReduced[reducedComponents[sourceVertex]] && reducedComponents[sourceVertex] == reducedComponents[targetVertex] && (state[boost::get(boost::edge_index, unreducedGraph, *currentEdge)] & OP_MASK))
			{
				Context::vertexPosition sourcePosition = vertexPositions[sourceVertex], targetPosition = vertexPositions[targetVertex];
				QGraphicsLineItem* newItem = new QGraphicsLineItem(sourcePosition.first, sourcePosition.second, targetPosition.first, targetPosition.second, reducedLinesItem);
				newItem->setPen(blackPen);
			}
		}
	}
}
