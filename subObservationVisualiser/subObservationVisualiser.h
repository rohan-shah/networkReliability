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
		subObservationVisualiser(NetworkReliabilitySubObsWithContext& subObsWithContextVector, float pointSize);
		~subObservationVisualiser();
		bool eventFilter(QObject* object, QEvent *event);
	private:
		void addBackgroundRectangle();
		void updateGraphics();
		void updateObservation();
		void fromStart();
		void addPoints();
		void addLines();

		void addReducedPoints();
		void addReducedLines();
		float pointSize;

		//Are we using a single sub-observation, or a vector of inputs?
		bool useSingleSubObs;
		QGraphicsScene* graphicsScene;
		QGraphicsView* graphicsView;
		QStatusBar* statusBar;
		QLabel* positionLabel;
		QLabel* reducedLabel;
		QFrame* statusFrame;
		QHBoxLayout* statusLayout;
		float minX, maxX, minY, maxY;
		//If we have a vector of input objects, which one are we currently looking at?
		int currentIndex;
		//and what is the current object? This is set for both single input / multiple inputs
		const NetworkReliabilitySubObs* subObs;
		//If the radius of the sub observation is 1, this holds the connected components of the sub observation. 
		std::vector<int> reducedComponents;
		//...and this holds the currently highlighted component
		int highlightedReducedComponent;
		//...and this hold whether or not we're looking at the reduced version or note. 
		bool reduced;
		//...and the reduced graph 
		NetworkReliabilitySubObs::getReducedGraphNoSelfWithWeightsInput reducedGraphData;
		//...and the number of components in the unreduced graph
		int nUnreducedComponents;
		const Context& context;
	};
}
#endif
