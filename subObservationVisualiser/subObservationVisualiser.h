#ifndef SUB_OBSERVATION_VISUALISER_HEADER_GUARD
#define SUB_OBSERVATION_VISUALISER_HEADER_GUARD
#include <QMainWindow>
#include "Context.h"
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QLabel>
#include <QStatusBar>
#include <QFrame>
#include <QHBoxLayout>
#include "NetworkReliabilityObs.h"
namespace networkReliability
{
	//If the next state is RESIMULATE, then we resimulate until we observe something that hits the next level, BUT
	class subObservationVisualiser : public QMainWindow
	{
		Q_OBJECT
	public:
		subObservationVisualiser(boost::shared_ptr<NetworkReliabilitySubObs> subObs, float pointSize);
		~subObservationVisualiser();
		bool eventFilter(QObject* object, QEvent *event);
	private:
		void addBackgroundRectangle();
		void updateGraphics();
		void fromStart();
		void addPoints();
		void addLines();

		void addReducedPoints();
		void addReducedLines();
		float pointSize;

		boost::shared_ptr<NetworkReliabilitySubObs> subObs;
		QGraphicsScene* graphicsScene;
		QGraphicsView* graphicsView;
		QStatusBar* statusBar;
		QLabel* positionLabel;
		QLabel* reducedLabel;
		QFrame* statusFrame;
		QHBoxLayout* statusLayout;
		float minX, maxX, minY, maxY;
		//If the radius of the sub observation is 1, this holds the connected components of the sub observation. 
		std::vector<int> radius1Components;
		//...and this holds the currently highlighted component
		int highlightedRadius1Component;
		//...and this hold whether or not we're looking at the reduced version or note. 
		bool reduced;
		//...and the reduced graph 
		NetworkReliabilitySubObs::getRadius1ReducedGraphNoSelfWithWeightsInput reducedGraphData;
		//...and the number of components in the unreduced graph
		int nUnreducedComponents;
	};
}
#endif
