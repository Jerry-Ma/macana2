#ifndef IMAGEVIEW_H
#define IMAGEVIEW_H

#include <QMainWindow>
#include <memory>

#include <nr3.h>
#include "colormapcontrol.h"
#include "qcustomplot/qcustomplot.h"

// macana-core types
class AnalParams;
class Array;
class Telescope;
class Observation;

namespace Ui {
class ImageView;
}

class ImageView : public QMainWindow
{
    Q_OBJECT

public:
    explicit ImageView(QWidget *parent = nullptr);
    ~ImageView();

    void initializeBeammaps(AnalParams* ap, Array* ar, Telescope* tel);

    static const std::vector<double>& getOrCreateColorMapDataCache(QCPColorMap*);

public slots:

    void showBeammap(int detId);
    void showStatus(QMouseEvent*);

private:
    Ui::ImageView *ui;

    // these are managed by mainwindow
    AnalParams* mAnalParams;
    Array* mArray;
    Telescope* mTelescope;
    // managed by imageview
    std::unique_ptr<Observation> mObservation;

    // beammap fit
    MatDoub fitParams;
    MatDoub previousFitParams;

    // plotting
    QMap<QString, QCPColorMap*> colorMaps;

    // plot controls
    QMap<QString, ColorMapControl*> colorMapControls;

    // panels
    QList<QString> colorMapNames {"signal", "model", "residual", "sigma", "weight"};

    void setColorMapData(QString name, double* data, std::function<double(double)> func);

    // cache for the plotted data
    static QMap<QCPColorMap*, std::vector<double>> colorMapDataCache;
};

#endif // IMAGEVIEW_H
