#include "subObservationVisualiserCollection.h"
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
	subObservationVisualiserCollection::subObservationVisualiserCollection(const NetworkReliabilitySubObsCollection& inputCollection, float pointSize)
		: collection(inputCollection), currentIndex(0)
	{
		if(collection.getSampleSize() == 0)
		{
			QTimer::singleShot(0, this, SLOT(close()));
		}
		else
		{
			statusBar = new subObservationStatusBar;
			setStatusBar(statusBar);

			base = new subObservationVisualiserBase(collection.getContext(), pointSize);
			base->installEventFilter(this);
			setCentralWidget(base);

			QObject::connect(base, &subObservationVisualiserBase::positionChanged, this, &subObservationVisualiserCollection::positionChanged);
			QObject::connect(base, &subObservationVisualiserBase::observationLeft, this, &subObservationVisualiserCollection::observationLeft);
			QObject::connect(base, &subObservationVisualiserBase::observationRight, this, &subObservationVisualiserCollection::observationRight);

			boost::shared_array<EdgeState> expandedState(new EdgeState[collection.getContext().getNEdges()]);
			collection.expand(currentIndex, expandedState);
			//Putting in dummy values for the last two constructor arguments
			NetworkReliabilitySubObs subObs(collection.getContext(), expandedState, collection.getRadius(), 0, 0);
			base->setObservation(subObs);
		}
	}
	subObservationVisualiserCollection::~subObservationVisualiserCollection()
	{
	}
	bool subObservationVisualiserCollection::eventFilter(QObject* object, QEvent *event)
	{
		if(event->type() == QEvent::KeyPress)
		{
			QKeyEvent* keyEvent = static_cast<QKeyEvent*>(event);
			if(keyEvent->key() == Qt::Key_R)
			{
				base->switchReduced();
				return true;
			}
			else if(keyEvent->key() == Qt::Key_Left && keyEvent->modifiers() & Qt::ShiftModifier)
			{
				observationLeft();
				return true;
			}
			else if(keyEvent->key() == Qt::Key_Right && keyEvent->modifiers() & Qt::ShiftModifier)
			{
				observationRight();
				return true;
			}
		}
		return false;
	}
	void subObservationVisualiserCollection::observationLeft()
	{
		if(currentIndex > 1)
		{
			currentIndex--;
			boost::shared_array<EdgeState> expandedState(new EdgeState[collection.getContext().getNEdges()]);
			collection.expand(currentIndex, expandedState);
			//Putting in dummy values for the last two constructor arguments
			NetworkReliabilitySubObs subObs(collection.getContext(), expandedState, collection.getRadius(), 0, 0);
			base->setObservation(subObs);
		}
	}
	void subObservationVisualiserCollection::observationRight()
	{
		if(currentIndex < (int)(collection.getSampleSize() - 1))
		{
			currentIndex++;
			boost::shared_array<EdgeState> expandedState(new EdgeState[collection.getContext().getNEdges()]);
			collection.expand(currentIndex, expandedState);
			//Putting in dummy values for the last two constructor arguments
			NetworkReliabilitySubObs subObs(collection.getContext(), expandedState, collection.getRadius(), 0, 0);
			base->setObservation(subObs);
		}
	}
	void subObservationVisualiserCollection::positionChanged(double x, double y)
	{
		statusBar->setPosition(x, y);
	}
}
