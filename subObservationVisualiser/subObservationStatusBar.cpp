#include "subObservationStatusBar.h"
#include <sstream>
namespace networkReliability
{
	subObservationStatusBar::subObservationStatusBar()
	{
		frame = new QFrame;
		addPermanentWidget(frame, 1);

		this->positionLabel = new QLabel;
		positionLabel->setText("");
		this->reducedLabel = new QLabel;
		reducedLabel->setText("");

		layout = new QHBoxLayout;
		layout->addWidget(positionLabel, 1, Qt::AlignLeft);
		layout->addWidget(reducedLabel, 0, Qt::AlignRight);
		layout->setContentsMargins(0,0,0,0);
		frame->setLayout(layout);
	}
	void subObservationStatusBar::setPosition(double x, double y)
	{
		std::stringstream ss;
		ss << "(" << x << ", " << y << ")";
		positionLabel->setText(QString::fromStdString(ss.str()));
	}
}
