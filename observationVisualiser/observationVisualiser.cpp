#include "observationVisualiser.h"
#include <QGraphicsView>
#include <QGraphicsItem>
#include <QEvent>
#include <QKeyEvent>
#include <boost/lexical_cast.hpp>
#include "graphAlgorithms.h"
#include "ZoomGraphicsView.h"
#include <QGraphicsSceneMouseEvent>
namespace networkReliability
{
	namespace observationVisualiserImpl
	{
		struct orderByFirst
		{
		public:
			bool operator()(::networkReliability::Context::vertexPosition const& first, ::networkReliability::Context::vertexPosition const& second) const
			{
				return first.first < second.first;
			}
		};
		struct orderBySecond
		{
			bool operator()(::networkReliability::Context::vertexPosition const& first, ::networkReliability::Context::vertexPosition const& second) const
			{
				return first.second < second.second;
			}
		};
	}
	observationVisualiser::observationVisualiser(Context const& context, boost::mt19937& randomSource, float pointSize)
		:obs(context, randomSource), context(context), randomSource(randomSource), pointSize(pointSize)
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

		this->statusLabel = new QLabel;
		this->positionLabel = new QLabel;
		statusLabel->setText("");
		positionLabel->setText("");

		statusLayout = new QHBoxLayout;
		statusLayout->addWidget(positionLabel, 1, Qt::AlignLeft);
		statusLayout->addWidget(statusLabel, 0, Qt::AlignRight);
		statusFrame->setLayout(statusLayout);
		setStatusBar(statusBar);

		const std::vector<Context::vertexPosition>& vertexPositions = context.getVertexPositions();
		observationVisualiserImpl::orderByFirst xOrder;
		observationVisualiserImpl::orderBySecond yOrder;
		minX = std::min_element(vertexPositions.begin(), vertexPositions.end(), xOrder)->first - pointSize;
		maxX = std::max_element(vertexPositions.begin(), vertexPositions.end(), xOrder)->first + pointSize;
		minY = std::min_element(vertexPositions.begin(), vertexPositions.end(), yOrder)->second - pointSize;
		maxY = std::max_element(vertexPositions.begin(), vertexPositions.end(), yOrder)->second + pointSize;

		setCentralWidget(graphicsView);
		updateGraphics();
	}
	void observationVisualiser::updateGraphics()
	{
		QList<QGraphicsItem*> allItems = graphicsScene->items();
		for(QList<QGraphicsItem*>::iterator i = allItems.begin(); i != allItems.end(); i++) delete *i;
		
		addBackgroundRectangle();
		addLines();
		addPoints();

		const EdgeState* state = obs.getState();

		std::vector<int> components;
		boost::detail::depth_first_visit_restricted_impl_helper<Context::internalGraph>::stackType stack;
		std::vector<boost::default_color_type> colorMap;
		countComponents(context, state, components, stack, colorMap);
		
		const std::vector<int>& interestVertices = context.getInterestVertices();
		
		std::vector<int> interestComponents;
		for(std::size_t i = 0; i < interestVertices.size(); i++)
		{
			interestComponents.push_back(components[interestVertices[i]]);
		}
		std::sort(interestComponents.begin(), interestComponents.end());
		interestComponents.erase(std::unique(interestComponents.begin(), interestComponents.end()), interestComponents.end());

		statusLabel->setText(QString::fromStdString("Interest vertices lie in " + boost::lexical_cast<std::string>(interestComponents.size()) + " components"));

		graphicsView->fitInView(minX, minY, maxX - minX, maxY - minY, Qt::KeepAspectRatio);
	}
	void observationVisualiser::addBackgroundRectangle()
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
	void observationVisualiser::addPoints()
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

		for(std::size_t vertexCounter = 0; vertexCounter < nVertices; vertexCounter++)
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
			QGraphicsSimpleTextItem* textItem = graphicsScene->addSimpleText(QString::fromStdString(boost::lexical_cast<std::string>(vertexCounter)));
			textItem->setBrush(redBrush);
			textItem->setPos(x, y);
		}
	}
	void observationVisualiser::addLines()
	{
		const EdgeState* state = obs.getState();
		const Context::internalGraph& graph = context.getGraph();
		const std::vector<Context::vertexPosition>& vertexPositions = context.getVertexPositions();

		QPen pen(QColor("black"));
		pen.setWidthF(pointSize/10);
		Context::internalGraph::edge_iterator start, end;
		boost::tie(start, end) = boost::edges(graph);

		while(start != end)
		{
			Context::vertexPosition sourcePosition = vertexPositions[start->m_source], targetPosition = vertexPositions[start->m_target];
			if(state[boost::get(boost::edge_index, graph, *start)] & OP_MASK)
			{
				graphicsScene->addLine(sourcePosition.first, sourcePosition.second, targetPosition.first, targetPosition.second, pen);
			}
			QGraphicsSimpleTextItem* text = graphicsScene->addSimpleText(QString::fromStdString(boost::lexical_cast<std::string>(boost::get(boost::edge_index, graph, *start))));
			text->setPos((sourcePosition.first + targetPosition.first)/2, (sourcePosition.second + targetPosition.second)/2);
			start++;
		}
	}
	bool observationVisualiser::eventFilter(QObject* object, QEvent *event)
	{
		if(event->type() == QEvent::GraphicsSceneMouseMove && object == graphicsScene)
		{
			QGraphicsSceneMouseEvent* mouseEvent = static_cast<QGraphicsSceneMouseEvent*>(event);
			QPointF position = mouseEvent->scenePos();
			std::stringstream ss;
			ss << "(" << position.x() << ", " << position.y() << ")";
			positionLabel->setText(QString::fromStdString(ss.str()));
		}
		else if(event->type() == QEvent::KeyPress)
		{
			QKeyEvent* keyEvent = static_cast<QKeyEvent*>(event);
			if(keyEvent->key() == Qt::Key_Enter || keyEvent->key() == Qt::Key_Return)
			{
				obs = NetworkReliabilityObs(context, randomSource);
				updateGraphics();
				return true;
			}
			return false;
		}
		return false;
	}
	observationVisualiser::~observationVisualiser()
	{
		delete statusLabel;
		delete positionLabel;
		delete statusLayout;

		delete statusFrame;
		delete statusBar;

		delete graphicsView;
		delete graphicsScene;
	}
}
