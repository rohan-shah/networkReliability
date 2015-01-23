#include "subObservationVisualiser.h"
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
	bool sortByFirst(Context::vertexPosition const& first, Context::vertexPosition const& second)
	{
		return first.first < second.first;
	}
	bool sortBySecond(Context::vertexPosition const& first, Context::vertexPosition const& second)
	{
		return first.second < second.second;
	}
	void subObservationVisualiser::updateReducedGraphData()
	{
		boost::detail::depth_first_visit_restricted_impl_helper<Context::internalGraph>::stackType stack;
		std::vector<boost::default_color_type> colorMap;
		if(useSingleSubObs)
		{
			nUnreducedComponents = countComponents(context, subObs->getState(), reducedComponents, stack, colorMap);
			subObs->getReducedGraphNoSelfWithWeights(reducedGraphData);
		}
		else
		{
			boost::shared_array<EdgeState> expandedState(new EdgeState[context.getNEdges()]);
			collection->expand(currentIndex, expandedState);
			//Entering dummy values for the last two inputs to the constructor
			NetworkReliabilitySubObs expandedSubObs(context, expandedState, collection->getRadius(), 0, 0);
			nUnreducedComponents = countComponents(context, expandedState.get(), reducedComponents, stack, colorMap);
			expandedSubObs.getReducedGraphNoSelfWithWeights(reducedGraphData);
		}
		updateGraphics();
	}
	void subObservationVisualiser::initialiseQt()
	{
		graphicsScene = new QGraphicsScene();
		graphicsScene->installEventFilter(this);
		graphicsScene->setItemIndexMethod(QGraphicsScene::NoIndex);

		graphicsView = new ZoomGraphicsView(graphicsScene);
		graphicsView->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
		graphicsView->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
		graphicsView->viewport()->installEventFilter(this);

		statusBar = new QStatusBar();
		statusFrame = new QFrame;
		statusBar->addPermanentWidget(statusFrame, 1);

		this->positionLabel = new QLabel;
		positionLabel->setText("");
		this->reducedLabel = new QLabel;
		reducedLabel->setText("");

		statusLayout = new QHBoxLayout;
		statusLayout->addWidget(positionLabel, 1, Qt::AlignLeft);
		statusLayout->addWidget(reducedLabel, 0, Qt::AlignRight);
		statusFrame->setLayout(statusLayout);

		setStatusBar(statusBar);

		const std::vector<Context::vertexPosition>& vertexPositions = context.getVertexPositions();
		minX = std::min_element(vertexPositions.begin(), vertexPositions.end(), sortByFirst)->first - pointSize;
		maxX = std::max_element(vertexPositions.begin(), vertexPositions.end(), sortByFirst)->first + pointSize;
		minY = std::min_element(vertexPositions.begin(), vertexPositions.end(), sortBySecond)->second - pointSize;
		maxY = std::max_element(vertexPositions.begin(), vertexPositions.end(), sortBySecond)->second + pointSize;

		setCentralWidget(graphicsView);
	}
	subObservationVisualiser::subObservationVisualiser(const NetworkReliabilitySubObsCollection& inputCollection, float pointSize)
		:pointSize(pointSize), useSingleSubObs(false), collection(&inputCollection), currentIndex(0), highlightedReducedComponent(-1), reduced(false), reducedGraphData(inputCollection.getContext().getInterestVertices()), context(inputCollection.getContext())
	{
		if(collection->getSampleSize() == 0)
		{
			QTimer::singleShot(0, this, SLOT(close()));
		}
		else
		{
			initialiseQt();
			observationChanged();
		}
	}
	void subObservationVisualiser::observationChanged()
	{
		highlightedReducedComponent = -1;
		updateReducedGraphData();
		updateGraphics();
	}
	subObservationVisualiser::subObservationVisualiser(const NetworkReliabilitySubObsWithContext& subObsWithContext, float pointSize)
		:pointSize(pointSize), useSingleSubObs(true), collection(NULL), currentIndex(0), subObs(&(subObsWithContext.getSubObs())), highlightedReducedComponent(-1), reduced(false), reducedGraphData(subObsWithContext.getSubObs().getContext().getInterestVertices()), context(subObsWithContext.getSubObs().getContext())
	{
		initialiseQt();
		observationChanged();
	}
	subObservationVisualiser::~subObservationVisualiser()
	{
	}
	void subObservationVisualiser::updateGraphics()
	{
		if(reduced)
		{
			reducedLabel->setText("Reduced");
		}
		else
		{
			reducedLabel->setText("Unreduced");
		}
		QList<QGraphicsItem*> allItems = graphicsScene->items();
		for(QList<QGraphicsItem*>::iterator i = allItems.begin(); i != allItems.end(); i++) delete *i;


		addBackgroundRectangle();
		if(!reduced)
		{
			addPoints();
			addLines();
		}
		else
		{
			addReducedPoints();
			addReducedLines();
		}
	}
	void subObservationVisualiser::addPoints()
	{
		std::size_t nVertices = boost::num_vertices(context.getGraph());

		const std::vector<Context::vertexPosition>& vertexPositions = context.getVertexPositions();
		const std::vector<int>& interestVertices = context.getInterestVertices();

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
			Context::vertexPosition currentPosition = vertexPositions[vertexCounter];
			float x = currentPosition.first;
			float y = currentPosition.second;
			if(highlightedReducedComponent != -1 && reducedComponents[vertexCounter] == highlightedReducedComponent)
			{
				graphicsScene->addEllipse(x - pointSize, y - pointSize, 2*pointSize, 2*pointSize, greenPen, greenBrush);
			}
			if(interestVertices.end() == std::find(interestVertices.begin(), interestVertices.end(), vertexCounter))
			{
				graphicsScene->addEllipse(x - pointSize/2, y - pointSize/2, pointSize, pointSize, blackPen, blackBrush);
			}
			else
			{
				graphicsScene->addEllipse(x - pointSize/2, y - pointSize/2, pointSize, pointSize, redPen, redBrush);
			}
		}
	}
	void subObservationVisualiser::addLines()
	{
		const EdgeState* state;
		boost::shared_array<EdgeState> expandedState;
		if(useSingleSubObs)
		{
			state = subObs->getState();
		}
		else
		{
			expandedState.reset(new EdgeState[context.getNEdges()]);
			collection->expand(currentIndex, expandedState);
			state = expandedState.get();
		}
		const Context::internalGraph& graph = context.getGraph();
		const std::vector<Context::vertexPosition>& vertexPositions = context.getVertexPositions();

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

		Context::internalGraph::edge_iterator start, end;
		boost::tie(start, end) = boost::edges(graph);

		boost::property_map<Context::internalGraph, boost::edge_index_t>::const_type edgeIndexMap = boost::get(boost::edge_index, graph);

		while(start != end)
		{
			Context::vertexPosition sourcePosition = vertexPositions[start->m_source], targetPosition = vertexPositions[start->m_target];
			if(state[boost::get(boost::edge_index, graph, *start)] & FIXED_OP)
			{
				graphicsScene->addLine(sourcePosition.first, sourcePosition.second, targetPosition.first, targetPosition.second, redPen);
			}
			else if(state[boost::get(boost::edge_index, graph, *start)] & UNFIXED_OP)
			{
				graphicsScene->addLine(sourcePosition.first, sourcePosition.second, targetPosition.first, targetPosition.second, pen);
			}
			else if(state[boost::get(boost::edge_index, graph, *start)] & UNFIXED_INOP)
			{
				graphicsScene->addLine(sourcePosition.first, sourcePosition.second, targetPosition.first, targetPosition.second, dashedPen);
			}
			else
			{
				graphicsScene->addLine(sourcePosition.first, sourcePosition.second, targetPosition.first, targetPosition.second, dashedRedPen);
			}
			QGraphicsSimpleTextItem* textItem = graphicsScene->addSimpleText(QString::fromStdString(boost::lexical_cast<std::string>(boost::get(edgeIndexMap, *start)))); 
			textItem->setPos((sourcePosition.first + targetPosition.first)/2, (sourcePosition.second + targetPosition.second)/2);
			start++;
		}
	}
	void subObservationVisualiser::addBackgroundRectangle()
	{
		QPen pen(Qt::NoPen);
		QColor grey("grey");
		grey.setAlphaF(0.5);

		QBrush brush;
		brush.setColor(grey);
		brush.setStyle(Qt::SolidPattern);
		QGraphicsRectItem* rect = graphicsScene->addRect(minX, minY, maxX - minX, maxY - minY, pen, brush);
		rect->setZValue(-1);
	}
	bool subObservationVisualiser::eventFilter(QObject* object, QEvent *event)
	{
		if(event->type() == QEvent::GraphicsSceneMouseMove && object == graphicsScene)
		{
			QGraphicsSceneMouseEvent* mouseEvent = static_cast<QGraphicsSceneMouseEvent*>(event);
			QPointF position = mouseEvent->scenePos();
			std::stringstream ss;
			ss << "(" << position.x() << ", " << position.y() << ")";
			positionLabel->setText(QString::fromStdString(ss.str()));
			const Context::internalGraph& graph = context.getGraph();
			std::size_t nVertices = boost::num_vertices(graph);
			const std::vector<Context::vertexPosition>& vertexPositions = context.getVertexPositions();

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
			else
			{
				
			}
		}
		else if(event->type() == QEvent::Leave && object == graphicsScene)
		{
			highlightedReducedComponent = -1;
			updateGraphics();
		}
		else if(event->type() == QEvent::KeyPress)
		{
			QKeyEvent* keyEvent = static_cast<QKeyEvent*>(event);
			if(keyEvent->key() == Qt::Key_R)
			{
				reduced = !reduced;
				if(!reduced) highlightedReducedComponent = -1;
				updateGraphics();
				return true;
			}
			else if(keyEvent->key() == Qt::Key_Left && !useSingleSubObs)
			{
				if(currentIndex > 1)
				{
					currentIndex--;
					observationChanged();
				}
			}
			else if(keyEvent->key() == Qt::Key_Right && !useSingleSubObs)
			{
				if(currentIndex < (int)(collection->getSampleSize() - 1))
				{
					currentIndex++;
					observationChanged();
				}
			}

		}
		return false;
	}
	void subObservationVisualiser::addReducedPoints()
	{
		const Context::internalGraph& unreducedGraph = context.getGraph();
		const std::vector<Context::vertexPosition>& vertexPositions = context.getVertexPositions();

		std::vector<bool> stillPresentInReduced(nUnreducedComponents, false);
		NetworkReliabilitySubObs::reducedGraphWithProbabilities::vertex_iterator currentReducedVertex, endReducedVertex;
		boost::tie(currentReducedVertex, endReducedVertex) = boost::vertices(reducedGraphData.outputGraph);
		for(;currentReducedVertex != endReducedVertex; currentReducedVertex++)
		{
			stillPresentInReduced[boost::get(boost::vertex_name, reducedGraphData.outputGraph, *currentReducedVertex)] = true;
		}
		
		QPen blackPen(QColor("black"));
		blackPen.setStyle(Qt::NoPen);
		QBrush blackBrush(QColor("black"));

		Context::internalGraph::vertex_iterator currentUnreducedVertex, endUnreducedVertex;
		boost::tie(currentUnreducedVertex, endUnreducedVertex) = boost::vertices(unreducedGraph);
		for(;currentUnreducedVertex != endUnreducedVertex; currentUnreducedVertex++)
		{
			if(stillPresentInReduced[reducedComponents[*currentUnreducedVertex]])
			{
				Context::vertexPosition currentPosition = vertexPositions[*currentUnreducedVertex];
				float x = currentPosition.first;
				float y = currentPosition.second;

				graphicsScene->addEllipse(x - pointSize/2, y - pointSize/2, pointSize, pointSize, blackPen, blackBrush);
			}
		}
	}
	void subObservationVisualiser::addReducedLines()
	{
		const Context::internalGraph& unreducedGraph = context.getGraph();
		const std::vector<Context::vertexPosition>& vertexPositions = context.getVertexPositions();

		std::vector<bool> stillPresentInReduced(nUnreducedComponents, false);
		NetworkReliabilitySubObs::reducedGraphWithProbabilities::vertex_iterator currentReducedVertex, endReducedVertex;
		boost::tie(currentReducedVertex, endReducedVertex) = boost::vertices(reducedGraphData.outputGraph);
		for(;currentReducedVertex != endReducedVertex; currentReducedVertex++)
		{
			stillPresentInReduced[boost::get(boost::vertex_name, reducedGraphData.outputGraph, *currentReducedVertex)] = true;
		}
		
		QPen blackPen(QColor("black"));

		Context::internalGraph::edge_iterator currentEdge, lastEdge;
		boost::tie(currentEdge, lastEdge) = boost::edges(unreducedGraph);
		const EdgeState* state;
		boost::shared_array<EdgeState> expandedState;
		if(useSingleSubObs) state = subObs->getState();
		else
		{
			expandedState.reset(new EdgeState[context.getNEdges()]);
			collection->expand(currentIndex, expandedState);
			state = expandedState.get();
		}
		for(; currentEdge != lastEdge; currentEdge++)
		{
			int sourceVertex = boost::source(*currentEdge, unreducedGraph);
			int targetVertex = boost::target(*currentEdge, unreducedGraph);
			if(stillPresentInReduced[reducedComponents[sourceVertex]] && reducedComponents[sourceVertex] == reducedComponents[targetVertex] && (state[boost::get(boost::edge_index, unreducedGraph, *currentEdge)] & OP_MASK))
			{
				Context::vertexPosition sourcePosition = vertexPositions[sourceVertex], targetPosition = vertexPositions[targetVertex];
				graphicsScene->addLine(sourcePosition.first, sourcePosition.second, targetPosition.first, targetPosition.second, blackPen);
			}
		}
	}
}
