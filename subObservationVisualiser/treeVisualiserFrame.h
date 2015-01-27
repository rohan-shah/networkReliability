#ifndef TREE_VISUALISER_GRAPHICS_VIEW_HEADER_GUARD
#define TREE_VISUALISER_GRAPHICS_VIEW_HEADER_GUARD
#include <QFrame>
#include <QGraphicsView>
#include <QGraphicsScene>
#include "NetworkReliabilitySubObsTree.h"
#include <QHBoxLayout>
#include <QGraphicsEllipseItem>
namespace networkReliability
{
	class treeVisualiserFrame : public QFrame
	{
		Q_OBJECT
	public:
		treeVisualiserFrame(const NetworkReliabilitySubObsTree& tree, float pointSize);
		void centreOn(double x, double y);
		void highlightPosition(double x, double y);
		bool eventFilter(QObject* object, QEvent *event);
	signals:
		void observationLeft();
		void observationRight();
		void observationUp();
		void observationDown();
		void vertexSelected(int vertex);
	private:
		double pointSize;
		QHBoxLayout* layout;
		QGraphicsView* graphicsView;
		QGraphicsScene* graphicsScene;
		QGraphicsEllipseItem* highlightItem;
	};
}
#endif
