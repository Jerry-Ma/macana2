#ifndef BEAMMAP_GUI_H
#define BEAMMAP_GUI_H

#include <QMainWindow>
#include "nr3.h"

class DomModel;  // model with parsed XML data

// macana-core types
class AnalParams;
class Array;
class TimePlace;
class Source;
class Telescope;
class Observation;

// plotting
class QCPColorMap;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

public slots:
    void openApXml();
    void openApXml(QString fp);

    void setObservation(int i);
    // void setObservation(QString f);

    void runBeammap();

    void showStatus(QMouseEvent*);
signals:
    void apLoaded(bool);

private:
    Ui::MainWindow *ui;

    // handles the input ap.xml file
    QString apXmlPath;
    DomModel* apXml;

    // macana-core objects
    std::unique_ptr<AnalParams> ap;
    int obsid;  // the current file index in the ap object
    std::unique_ptr<Array> array;
    std::unique_ptr<TimePlace> timePlace;
    std::unique_ptr<Source> source;
    std::unique_ptr<Telescope> telescope;
    std::unique_ptr<Observation> observation;

    // beammap data
    MatDoub fitParams;
    MatDoub previousFitParams;
    VecDoub needsIteration;
    MatDoub originalhVals;
    int iteration = 0;
    int cap = 0;
    double cutOff = 0.;

    // plotting
    QCPColorMap* colorMap;
    int did = 0;  // detector index to plot
    void showMap(int did);
};

#endif // BEAMMAP_GUI_H
