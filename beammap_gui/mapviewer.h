#ifndef MAPVIEWER_H
#define MAPVIEWER_H

#include <QWidget>

namespace Ui {
class MapViewer;
}

class MapViewer : public QWidget
{
    Q_OBJECT

public:
    explicit MapViewer(QWidget *parent = nullptr);
    ~MapViewer();

private:
    Ui::MapViewer *ui;
};

#endif // MAPVIEWER_H
