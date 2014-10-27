#ifndef SPLITTING_VISUALISER_HEADER_GUARD
#define SPLITTING_VISUALISER_HEADER_GUARD
#include <QMainWindow>
#include "Context.h"
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QLabel>
#include <QStatusBar>
#include "NetworkReliabilityObs.h"
namespace networkReliability
{
	enum nextState
	{
		DECREASE_RADIUS, RESIMULATE
	};
	//If the next state is RESIMULATE, then we resimulate until we observe something that hits the next level, BUT
	class splittingVisualiser : public QMainWindow
	{
		Q_OBJECT
	public:
		splittingVisualiser(Context const& context, int seed, float pointSize, int initialRadius);
		~splittingVisualiser();
		bool eventFilter(QObject* object, QEvent *event);
	private:
		void addBackgroundRectangle();
		void updateGraphics(int connectionRadius, int highlightRadius);
		void fromStart();
		void nextStep();
		void addPoints();
		void addLines();
		Context const& context;
		boost::mt19937 randomSource;
		float pointSize;
		int seed;
		int initialRadius;
		int currentRadius;

		NetworkReliabilityObs obs;
		QGraphicsScene* graphicsScene;
		QGraphicsView* graphicsView;
		QStatusBar* statusBar;
		QLabel* statusLabel;

		nextState nextAction;
		float minX, maxX, minY, maxY;
	};
}
#endif