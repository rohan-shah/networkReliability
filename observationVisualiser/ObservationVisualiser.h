#ifndef NETWORK_RELIABILITY_HEADER_GUARD
#define NETWORK_RELIABILITY_HEADER_GUARD
#include <QMainWindow>
#include "NetworkReliabilityObs.h"
#include <QGraphicsView>
#include <QGraphicsScene>
#include <QStatusBar>
#include <QLabel>
#include <boost/random/mersenne_twister.hpp>
#include <boost/shared_ptr.hpp>
namespace networkReliability
{
	class ObservationVisualiser : public QMainWindow
	{
		Q_OBJECT
	public:
		ObservationVisualiser(Context const& context, boost::mt19937& randomSource, float pointSize);
		~ObservationVisualiser();
		bool eventFilter(QObject* object, QEvent *event);
	private:
		//the different graphics bits that have to get added
		void addBackgroundRectangle();
		void addPoints();
		void addLines();
		void updateGraphics();

		NetworkReliabilityObs obs;
		Context const& context;
		boost::mt19937& randomSource;
		QGraphicsScene* graphicsScene;
		QGraphicsView* graphicsView;
		QStatusBar* statusBar;
		QLabel* statusLabel;
		float minX, maxX, minY, maxY;
		float pointSize;
		float probability;
	};
}
#endif