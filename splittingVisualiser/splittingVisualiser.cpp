#include "splittingVisualiser.h"
#include <QEvent>
#include <QKeyEvent>
#include "ZoomGraphicsView.h"
#include <QGraphicsRectItem>
#include "graphAlgorithms.h"
#include <boost/lexical_cast.hpp>
namespace networkReliability
{
	namespace splittingVisualiserImpl
	{
		template <typename Proj> struct less_by_t 
		{
			Proj p;
			template <typename T>
			bool operator ()(T const& a, T const& b) const {
				return p(a) < p(b);
			};
		};
		template <typename Proj> less_by_t<Proj> less_by(Proj p) 
		{
			less_by_t<Proj> ret;
			ret.p = p;
			return ret;
		};
	}
	splittingVisualiser::splittingVisualiser(Context const& context, int seed, float pointSize, int initialRadius)
		:context(context), pointSize(pointSize), initialRadius(initialRadius), obs(context, randomSource), currentRadius(initialRadius), seed(seed)
	{
		randomSource.seed(seed);
		if(initialRadius < 0)
		{
			throw std::runtime_error("Internal error");
		}
		graphicsScene = new QGraphicsScene();
		graphicsScene->installEventFilter(this);
		graphicsScene->setItemIndexMethod(QGraphicsScene::NoIndex);

		graphicsView = new ZoomGraphicsView(graphicsScene);
		graphicsView->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
		graphicsView->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
		graphicsView->viewport()->installEventFilter(this);
		
		statusBar = new QStatusBar();
		this->statusLabel = new QLabel;
		statusLabel->setText("");
		statusBar->addPermanentWidget(statusLabel);
		setStatusBar(statusBar);
		const std::vector<Context::vertexPosition>& vertexPositions = context.getVertexPositions();


		auto xLambda = [](Context::vertexPosition const& x){ return x.first;};
		auto yLambda = [](Context::vertexPosition const& x){ return x.second;};
		auto xPredicate = splittingVisualiserImpl::less_by(xLambda);
		auto yPredicate = splittingVisualiserImpl::less_by(yLambda);
		minX = std::min_element(vertexPositions.begin(), vertexPositions.end(), xPredicate)->first - pointSize;
		maxX = std::max_element(vertexPositions.begin(), vertexPositions.end(), xPredicate)->first + pointSize;
		minY = std::min_element(vertexPositions.begin(), vertexPositions.end(), yPredicate)->second - pointSize;
		maxY = std::max_element(vertexPositions.begin(), vertexPositions.end(), yPredicate)->second + pointSize;

		setCentralWidget(graphicsView);

		fromStart();
	}
	void splittingVisualiser::fromStart()
	{
		seed++;
		randomSource.seed(seed);

		const std::size_t nEdges = boost::num_edges(context.getGraph());
		std::vector<int> components(nEdges);
		std::vector<boost::default_color_type> colorMap;
		boost::detail::depth_first_visit_restricted_impl_helper<Context::internalGraph>::stackType stack;
		currentRadius = initialRadius;
		while(true)
		{
			obs = NetworkReliabilityObs::constructConditional(context, randomSource);
			NetworkReliabilitySubObs subObs = obs.getSubObservation(initialRadius);
			if(!isSingleComponent(context, subObs.getState(), components, stack, colorMap)) break;
		}
		updateGraphics();
	}
	splittingVisualiser::~splittingVisualiser()
	{
	}
	void splittingVisualiser::updateGraphics()
	{
		QList<QGraphicsItem*> allItems = graphicsScene->items();
		for(QList<QGraphicsItem*>::iterator i = allItems.begin(); i != allItems.end(); i++) delete *i;

		addBackgroundRectangle();
		addPoints();
		addLines();
		statusLabel->setText(QString::fromStdString("Current observation is potentially disconnected at radius " + boost::lexical_cast<std::string>(currentRadius)));
	}
	void splittingVisualiser::addPoints()
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

		for(int vertexCounter = 0; vertexCounter < nVertices; vertexCounter++)
		{
			Context::vertexPosition currentPosition = vertexPositions[vertexCounter];
			float x = currentPosition.first;
			float y = currentPosition.second;
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
	void splittingVisualiser::addLines()
	{
		const EdgeState* state = obs.getState();
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
	void splittingVisualiser::addBackgroundRectangle()
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
	void splittingVisualiser::nextStep()
	{
		const Context::internalGraph& graph = context.getGraph();
		const std::size_t nEdges = boost::num_edges(graph);

		typedef boost::color_traits<boost::default_color_type> Color;
		//define the stuff used in the connected components algorithm just once, and reuse it. 
		std::vector<int> components(nEdges);
		std::vector<boost::default_color_type> colorMap;
		boost::detail::depth_first_visit_restricted_impl_helper<Context::internalGraph>::stackType stack;
		if(currentRadius > 0)
		{
			NetworkReliabilitySubObs subObs = obs.getSubObservation(currentRadius);
			if(subObs.getMinCut() >= HIGH_CAPACITY)
			{
				throw std::runtime_error("Mincut was infeasibly large");
			}
			currentRadius--;
			while(true)
			{
				obs = subObs.getObservation(randomSource);
				NetworkReliabilitySubObs newSubObs = obs.getSubObservation(currentRadius);
				if(isSingleComponent(context, newSubObs.getState(), components, stack, colorMap)) continue;
				updateGraphics();
				if(newSubObs.getMinCut() >= HIGH_CAPACITY)
				{
					throw std::runtime_error("Mincut was infeasibly large");
				}
				break;
			}
		}
	}
	bool splittingVisualiser::eventFilter(QObject*, QEvent *event)
	{
		if(event->type() == QEvent::KeyPress)
		{
			QKeyEvent* keyEvent = static_cast<QKeyEvent*>(event);
			if(keyEvent->key() == Qt::Key_Enter || keyEvent->key() == Qt::Key_Return)
			{
				nextStep();
			}
			if(keyEvent->key() == Qt::Key_N)
			{
				fromStart();
			}
		}
		return false;
	}
}