#include "subObservationVisualiserCollection.h"
#include <QEvent>
#include <QKeyEvent>
#include <QTimer>
#include "zoomGraphicsView.h"
#include <QGraphicsRectItem>
#include "graphAlgorithms.h"
#include <boost/lexical_cast.hpp>
#include <QGraphicsSceneMouseEvent>
namespace networkReliability
{
	subObservationVisualiserCollection::subObservationVisualiserCollection(const NetworkReliabilityObsCollection& inputCollection, float pointSize, std::string text)
		: collection(inputCollection), currentIndex(0), text(text)
	{
		if(collection.getSampleSize() == 0)
		{
			QTimer::singleShot(0, this, SLOT(close()));
		}
		else
		{
			statusBar = new subObservationStatusBar(text);
			setStatusBar(statusBar);

			base = new subObservationVisualiserBase(collection.getContext(), pointSize);
			base->installEventFilter(this);
			setCentralWidget(base);

			QObject::connect(base, &subObservationVisualiserBase::positionChanged, this, &subObservationVisualiserCollection::positionChanged);
			QObject::connect(base, &subObservationVisualiserBase::observationLeft, this, &subObservationVisualiserCollection::observationLeft);
			QObject::connect(base, &subObservationVisualiserBase::observationRight, this, &subObservationVisualiserCollection::observationRight);

			boost::shared_array<edgeState> expandedState(new edgeState[collection.getContext().getNEdges()]);
			collection.expand(currentIndex, expandedState);
			//Putting in dummy values for the last two constructor arguments
			NetworkReliabilityObs subObs(collection.getContext(), expandedState);
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
		if(currentIndex > 0)
		{
			currentIndex--;
			boost::shared_array<edgeState> expandedState(new edgeState[collection.getContext().getNEdges()]);
			collection.expand(currentIndex, expandedState);
			//Putting in dummy values for the last two constructor arguments
			NetworkReliabilityObs obs(collection.getContext(), expandedState);
			base->setObservation(obs);
		}
	}
	void subObservationVisualiserCollection::observationRight()
	{
		if(currentIndex < (int)(collection.getSampleSize() - 1))
		{
			currentIndex++;
			boost::shared_array<edgeState> expandedState(new edgeState[collection.getContext().getNEdges()]);
			collection.expand(currentIndex, expandedState);
			//Putting in dummy values for the last two constructor arguments
			NetworkReliabilityObs obs(collection.getContext(), expandedState);
			base->setObservation(obs);
		}
	}
	void subObservationVisualiserCollection::positionChanged(double x, double y)
	{
		statusBar->setPosition(x, y);
	}
}
