#include "subObservationVisualiserSingle.h"
#include <QEvent>
#include <QKeyEvent>
#include <QTimer>
#include "ZoomGraphicsView.h"
#include <QGraphicsRectItem>
#include "graphAlgorithms.h"
#include <boost/lexical_cast.hpp>
#include <QGraphicsSceneMouseEvent>
namespace networkReliability
{
	subObservationVisualiserSingle::subObservationVisualiserSingle(const NetworkReliabilitySubObsWithContext& subObsWithContext, float pointSize)
		:subObsWithContext(subObsWithContext) 
	{
		statusBar = new subObservationStatusBar;
		setStatusBar(statusBar);
		
		base = new subObservationVisualiserBase(subObsWithContext.getContext(), pointSize);
		base->installEventFilter(this);

		setCentralWidget(base);
		QObject::connect(base, &subObservationVisualiserBase::positionChanged, this, &subObservationVisualiserSingle::positionChanged);
		base->setObservation(subObsWithContext.getSubObs());
	}
	void subObservationVisualiserSingle::positionChanged(double x, double y)
	{
		statusBar->setPosition(x, y);
	}
	subObservationVisualiserSingle::~subObservationVisualiserSingle()
	{}
	bool subObservationVisualiserSingle::eventFilter(QObject* object, QEvent *event)
	{
		if(event->type() == QEvent::KeyPress)
		{
			QKeyEvent* keyEvent = static_cast<QKeyEvent*>(event);
			if(keyEvent->key() == Qt::Key_R)
			{
				base->switchReduced();
				return true;
			}
		}
		return false;
	}
}