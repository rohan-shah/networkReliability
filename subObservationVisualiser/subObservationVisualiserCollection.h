#ifndef SUB_OBSERVATION_VISUALISER_COLLECTION_HEADER_GUARD
#define SUB_OBSERVATION_VISUALISER_COLLECTION_HEADER_GUARD
#include <QMainWindow>
#include "Context.h"
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QLabel>
#include <QStatusBar>
#include <QFrame>
#include <QHBoxLayout>
#include "subObservationVisualiserBase.h"
#include "subObservationStatusBar.h"
#include "NetworkReliabilitySubObsCollection.h"
namespace networkReliability
{
	//If the next state is RESIMULATE, then we resimulate until we observe something that hits the next level, BUT
	class subObservationVisualiserCollection : public QMainWindow
	{
		Q_OBJECT
	public:
		subObservationVisualiserCollection(const NetworkReliabilitySubObsCollection& collection, float pointSize, std::string text = "");
		~subObservationVisualiserCollection();
		bool eventFilter(QObject* object, QEvent *event);
	public slots:
		void positionChanged(double x, double y);
		void observationLeft();
		void observationRight();
	private:
		subObservationStatusBar* statusBar;
		subObservationVisualiserBase* base;
		const NetworkReliabilitySubObsCollection& collection;
		int currentIndex;
		std::string text;
	};
}
#endif
