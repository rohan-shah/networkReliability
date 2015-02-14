#ifndef SUB_OBSERVATION_VISUALISER_SINGLE_HEADER_GUARD
#define SUB_OBSERVATION_VISUALISER_SINGLE_HEADER_GUARD
#include <QMainWindow>
#include "Context.h"
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QLabel>
#include <QStatusBar>
#include <QFrame>
#include <QHBoxLayout>
#include "NetworkReliabilityObs.h"
#include "subObservationVisualiserBase.h"
#include "subObservationStatusBar.h"
namespace networkReliability
{
	class subObservationVisualiserSingle : public QMainWindow
	{
		Q_OBJECT
	public:
		subObservationVisualiserSingle(const NetworkReliabilityObsWithContext& obsWithContext, float pointSize);
		~subObservationVisualiserSingle();
		bool eventFilter(QObject* object, QEvent *event);
	public slots:
		void positionChanged(double x, double y);
	private:
		subObservationStatusBar* statusBar;
		subObservationVisualiserBase* base;
		const NetworkReliabilityObsWithContext& obsWithContext;
	};
}
#endif
