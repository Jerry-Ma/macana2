#include "beammap_gui/mapviewer.h"
#include "ui_mapviewer.h"

MapViewer::MapViewer(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::MapViewer)
{
    ui->setupUi(this);
}

MapViewer::~MapViewer()
{
    delete ui;
}
