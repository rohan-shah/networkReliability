#ifndef SUB_OBSERVATION_VISUALISER_HEADER_GUARD
#define SUB_OBSERVATION_VISUALISER_HEADER_GUARD
#include <QMainWindow>
#include "Context.h"
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QLabel>
#include <QStatusBar>
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
		float pointSize;

		boost::shared_ptr<NetworkReliabilitySubObs> subObs;
		QGraphicsScene* graphicsScene;
		QGraphicsView* graphicsView;

		float minX, maxX, minY, maxY;
	};
}
#endif