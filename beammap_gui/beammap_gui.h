#ifndef BEAMMAP_GUI_H
#define BEAMMAP_GUI_H

#include <QMainWindow>
#include <nr3.h>

class DomModel;  // model with parsed XML data

// macana-core types
class AnalParams;
class Array;
class TimePlace;
class Source;
class Telescope;
class Observation;

// class MatDoub;
// class VecDoub;

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
    void openApXml(const QString& path);
    void setObsId(int idx);

    void runBeammap();

    void showStatus(QMouseEvent*);
signals:
    void observationInitialized(bool);

private:
    Ui::MainWindow *ui;

    // handles the input ap.xml file
    QString apXmlPath;
    DomModel* apXml = nullptr;
    int obsId = 0;  // the current file index in the ap.xml

    void initializeObservation();  // initialize macana-core objects with optional file number set
    // macana-core
    std::unique_ptr<AnalParams> mAnalParams;
    std::unique_ptr<Array> mArray;
    std::unique_ptr<TimePlace> mTimePlace;
    std::unique_ptr<Source> mSource;
    std::unique_ptr<Telescope> mTelescope;
    std::unique_ptr<Observation> mObservation;

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
