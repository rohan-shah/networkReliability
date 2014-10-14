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
	class splittingVisualiser : public QMainWindow
	{
		Q_OBJECT
	public:
		splittingVisualiser(Context const& context, int seed, float pointSize, int initialRadius);
		~splittingVisualiser();
		bool eventFilter(QObject* object, QEvent *event);
	private:
		void addBackgroundRectangle();
		void updateGraphics();
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

		float minX, maxX, minY, maxY;
	};
}
#endif