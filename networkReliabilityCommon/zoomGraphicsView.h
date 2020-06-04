#ifndef ZOOM_GRAPHICS_HEADER_GUARD
#define ZOOM_GRAPHICS_HEADER_GUARD
#include <QGraphicsView>
namespace networkReliability
{
	class zoomGraphicsView : public QGraphicsView
	{
	public:
		zoomGraphicsView(QGraphicsScene* scene);
	protected:
		void wheelEvent(QWheelEvent* event);
		void keyPressEvent(QKeyEvent* event);
	};
}
#endif
