#ifndef SUB_OBSERVATION_VISUALISER_BASE_HEADER_GUARD
#define SUB_OBSERVATION_VISUALISER_BASE_HEADER_GUARD
#include <QFrame>
#include "Context.h"
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QLabel>
#include "NetworkReliabilitySubObs.h"
#include <QHBoxLayout>
namespace networkReliability
{
	class subObservationVisualiserBase : public QFrame
	{
		Q_OBJECT
	public:
		subObservationVisualiserBase(const Context& context, float pointSize);
		~subObservationVisualiserBase();
		bool eventFilter(QObject* object, QEvent *event);
		void setReduced(bool reduced);
		void switchReduced();
		void setObservation(const NetworkReliabilitySubObs& subObs);
	signals:
		void positionChanged(double x, double y);
		void observationLeft();
		void observationRight();
		void observationUp();
		void observationDown();
	private:
		void initialiseControls();
		void addBackgroundRectangle();
		void constructGraphics(const NetworkReliabilitySubObs& subObs);
		void updateGraphics();
		void updateReducedGraphData(const NetworkReliabilitySubObs& subObs);
		void constructUnreducedPoints();
		void constructUnreducedLines(const NetworkReliabilitySubObs& subObs);

		void constructReducedPoints();
		void constructReducedLines(const NetworkReliabilitySubObs& subObs);
		float pointSize;

		QGraphicsScene* graphicsScene;
		QGraphicsView* graphicsView;
		QHBoxLayout* layout;
		float minX, maxX, minY, maxY;
		//If the radius of the sub observation is 1, this holds the connected components of the sub observation. 
		std::vector<int> reducedComponents;
		//...and this holds the currently highlighted component
		int highlightedReducedComponent;
		//...and the reduced graph 
		NetworkReliabilitySubObs::getReducedGraphNoSelfWithWeightsInput reducedGraphData;
		//...and the number of components in the unreduced graph
		int nUnreducedComponents;
		bool reduced;
		const Context& context;

		QGraphicsItemGroup* reducedPointsItem;
		QGraphicsItemGroup* reducedLinesItem;
		QGraphicsItemGroup* unreducedPointsItem;
		QGraphicsItemGroup* unreducedLinesItem;
		QGraphicsItemGroup* highlightedGroupItem;
		QGraphicsRectItem* backgroundItem;
	};
}
#endif
