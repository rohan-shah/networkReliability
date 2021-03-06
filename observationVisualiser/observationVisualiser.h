#ifndef NETWORK_RELIABILITY_HEADER_GUARD
#define NETWORK_RELIABILITY_HEADER_GUARD
#include <QMainWindow>
#include "networkReliabilityObs.h"
#include <QGraphicsView>
#include <QGraphicsScene>
#include <QStatusBar>
#include <QLabel>
#include <QFrame>
#include <QHBoxLayout>
#include <boost/random/mersenne_twister.hpp>
#include <boost/shared_ptr.hpp>
namespace networkReliability
{
	class observationVisualiser : public QMainWindow
	{
		Q_OBJECT
	public:
		observationVisualiser(context const& contextObj, boost::mt19937& randomSource, float pointSize);
		~observationVisualiser();
		bool eventFilter(QObject* object, QEvent *event);
	private:
		//the different graphics bits that have to get added
		void addBackgroundRectangle();
		void addPoints();
		void addLines();
		void updateGraphics();

		NetworkReliabilityObs obs;
		context const& contextObj;
		boost::mt19937& randomSource;
		QGraphicsScene* graphicsScene;
		QGraphicsView* graphicsView;
		QStatusBar* statusBar;
		QFrame* statusFrame;
		QHBoxLayout* statusLayout;
		QLabel* positionLabel;
		QLabel* statusLabel;
		float minX, maxX, minY, maxY;
		float pointSize;
		float probability;
	};
}
#endif
