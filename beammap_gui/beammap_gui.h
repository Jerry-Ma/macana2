#ifndef BEAMMAP_GUI_H
#define BEAMMAP_GUI_H

#include <QMainWindow>
#include <QStringListModel>

#include <nr3.h>
#include "dommodel.h"

// macana-core types
class AnalParams;
class Array;
class TimePlace;
class Source;
class Telescope;
class Observation;

// class MatDoub;
// class VecDoub;


namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow() override;
    void showEvent(QShowEvent* event) override;

public slots:
    void openApXml();
    void openApXml(const QString& path);
    void setObsId(int idx);
    void setDetId(int idx);

    void runBeammap();

signals:
    void observationInitialized();

private:
    Ui::MainWindow *ui;
    bool startUp = true;

    // handles the input ap.xml file
    QString apXmlPath;
    DomModel* apXml = nullptr;
    QStringListModel detList;
    int obsId = 0;  // the current file index in the ap.xml
    int detId = 0;  // detector index to work with
    void updateDetectorIndices();  // update di in mArray as well as the detComboBox

    void initializeObservation();  // initialize macana-core objects with optional file number set
    // macana-core
    std::unique_ptr<AnalParams> mAnalParams;
    std::unique_ptr<Array> mArray;
    std::unique_ptr<TimePlace> mTimePlace;
    std::unique_ptr<Source> mSource;
    std::unique_ptr<Telescope> mTelescope;
    std::unique_ptr<Observation> mObservation;

    // beammap data
    VecDoub needsIteration;
    MatDoub originalhVals;

    int iteration = 0;
    int cap = 0;
    double cutOff = 0.;


};

#endif // BEAMMAP_GUI_H
