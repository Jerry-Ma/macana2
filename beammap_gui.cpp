#include <QApplication>
#include <QDomDocument>
#include "beammap_gui.h"
#include "beammap_gui/dommodel.h"
#include "ui_beammap_gui.h"

#include "AnalParams.h"
#include "Observation.h"
#include "TimePlace.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::openApXml()
{
    QString filePath = QFileDialog::getOpenFileName(this, tr("Open File"),
        xmlPath, tr("XML files (*.xml)"));

    if (!filePath.isEmpty()) {
        QFile file(filePath);
        if (file.open(QIODevice::ReadOnly)) {
            QDomDocument document;
            if (document.setContent(&file)) {
                DomModel *newModel = new DomModel(document, this);
                ui->treeView->setModel(newModel);
                delete model;
                model = newModel;
                xmlPath = filePath;
            }
            file.close();
        }
    }
}

void MainWindow::runBeammap()
{
    if (xmlPath == "") {
        openApXml();
    }
    AnalParams* ap = new AnalParams(xmlPath.toStdString(), 1);
    int nFiles = ap->getNFiles();
    qDebug() << "number of ap files: " << nFiles;
    //begin loop over input files
    for(int fileNum = 0; fileNum < nFiles; ++fileNum){
        qDebug() << "starting file number " << fileNum;
        ap->setDataFile(fileNum);
        TimePlace* timePlace = new TimePlace(ap);
        Array* array = new Array(ap);
        array->populate();
        Source* source = new Source(ap, timePlace);
        Telescope* telescope = new Telescope(ap, timePlace, source);
        array->updateDetectorIndices();
        int *di = array->getDetectorIndices();

        //setting up physical coordinate system
        telescope->absToPhysEqPointing();
        for(int i = 0; i < array->getNDetectors(); i++){
          array->detectors[di[i]].getPointing(telescope, timePlace, source);
          array->detectors[di[i]].getAzElPointing(telescope);
        }
        array->findMinMaxXY();

        //despiking at detector level, detector by detector
        for(int i = 0; i < array->getNDetectors(); i++){
          array->detectors[di[i]].despike(ap->getDespikeSigma());
        }
        array->updateDetectorIndices();
        di = array->getDetectorIndices();

        //replace flagged data in scans with faked data (useful for pca if needed)
        array->fakeFlaggedData(telescope);

        //make a fake source for psf determination
        for(int i = 0; i <array->getNDetectors(); i++){
          array->detectors[di[i]].makeKernelTimestream(telescope);
        }

        //lowpass the data
        for(int i = 0; i < array->getNDetectors(); i++){
          array->detectors[di[i]].lowpass(&array->digFiltTerms[0], array->nFiltTerms);
        }

        VecBool obsFlags(array->detectors[0].getNSamples());
        for(int i = 0; i < array->detectors[0].getNSamples(); i++){
          obsFlags[i] = 0;
          for(int j = 0; j < array->getNDetectors(); j++){
            if(array->detectors[di[j]].hSampleFlags[i]) obsFlags[i] = 1;
          }
        }
        telescope->checkScanLengths(obsFlags, array->detectors[0].getSamplerate());

        //generate the pointing signals
        for(int i=0;i<array->getNDetectors();i++){
          array->detectors[di[i]].estimateExtinction(array->getAvgTau());
        }

        //create vectors to store future fit parameters for iterative cleaning
        MatDoub fitParams(array->getNDetectors(), 7, 0.);
        MatDoub previousFitParams(array->getNDetectors(), 7, 1.0);

        //create the vector that stores whether to terminate iteration for each detector
        VecDoub needsIteration(array->getNDetectors(), 1.0);

        //store original hValues for iterative cleaning
        int nSamples = array->detectors[di[0]].getNSamples();
        MatDoub originalhVals(array->getNDetectors(), nSamples, 0.);
        for(int i = 0; i < array->getNDetectors(); i++){
          for(int j = 0; j < nSamples; j++){
            originalhVals[i][j] = array->detectors[di[i]].hValues[j];
          }
        }
    }
}

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();

    return a.exec();
}