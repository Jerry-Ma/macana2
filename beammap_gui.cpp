#include <QApplication>
#include <QDomDocument>
#include "beammap_gui.h"
#include "beammap_gui/dommodel.h"
#include "beammap_gui/domitem.h"
#include "ui_beammap_gui.h"

#include "AnalParams.h"
#include "Observation.h"
#include "TimePlace.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    colorMap = new QCPColorMap(ui->canvas->xAxis, ui->canvas->yAxis);
    colorMap->setInterpolate(false);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::openApXml()
{
     QString fp = QFileDialog::getOpenFileName(
                this, tr("Open File"),
                apXmlPath, tr("XML files (*.xml)"));
     return openApXml(fp);
}

void MainWindow::openApXml(QString fp)
{
    qDebug() << "load" << fp;
    apXmlPath = fp;
    // initialize the macana-core objects
    ap = std::make_unique<AnalParams>(apXmlPath.toStdString(), 1);
    qDebug() << "number of observations:" << ap->getNFiles();
    if (!fp.isEmpty()) {
        QFile fo(fp);
        if (fo.open(QIODevice::ReadWrite)) {
            QDomDocument doc;
            if (doc.setContent(&fo)) {
                // setup ap.xml view
                auto model = new DomModel(doc, this);
                ui->apXmlTreeView->setModel(model);
                ui->apXmlTreeView->expandAll();
                ui->obsTreeView->setModel(model);
                ui->obsTreeView->setRootIndex(QModelIndex());
                ui->obsTreeView->expandAll();
                ui->obsComboBox->blockSignals(true); // silence until we finish setup
                // set up a proxy to populate the observation combobox
                QSortFilterProxyModel* proxyModel = new QSortFilterProxyModel(this);
                proxyModel->setSourceModel(model);
                // keep only the file entries
                proxyModel->setRecursiveFilteringEnabled(true);
                proxyModel->setFilterRegExp(QRegExp(
                            "f\\d+",
                            Qt::CaseInsensitive));
                ui->obsComboBox->setModel(proxyModel);
                // get the index for observation tag
                QModelIndex obsindex = proxyModel->match(
                        proxyModel->index(0, 0),
                        Qt::DisplayRole,
                        QVariant::fromValue(QString("observations")),
                        1,
                        Qt::MatchExactly | Qt::MatchRecursive
                        ).at(0);
                ui->obsComboBox->setRootModelIndex(obsindex);
                ui->obsComboBox->blockSignals(false); // let go
                ui->obsComboBox->setCurrentIndex(0);
                // delete apXml;  // cleanup the old model
                apXml = model;
            }
            fo.close();
        }
    }
    emit apLoaded(true);
}

void MainWindow::setObservation(int obsid)
{
    this->obsid = obsid;
    // figure out the model index
    auto proxyModel = static_cast<QSortFilterProxyModel*>(ui->obsComboBox->model());
    QModelIndex obsindex = proxyModel->match(
        proxyModel->index(0, 0),
        Qt::DisplayRole,
        QVariant::fromValue(QString("observations")),
        1,
        Qt::MatchExactly | Qt::MatchRecursive
        ).at(0);
    QModelIndex currentindex = proxyModel->mapToSource(
                    proxyModel->index(obsid, 0, obsindex));
    // set the obs tree view to current obs
    ui->obsTreeView->setRootIndex(currentindex);
    DomItem* obsitem = static_cast<DomItem*>(currentindex.internalPointer());
    qDebug() << "work on file #" << obsid
             << ui->obsComboBox->currentText()
             << obsitem->node().firstChildElement().text();

    // setup macana-core objects
    ap->setDataFile(obsid);
    // create an Array object
    array = std::make_unique<Array>(ap.get());
    array->populate();
    // make a TimePlace
    timePlace = std::make_unique<TimePlace>(ap.get());
    // make a Source
    source = std::make_unique<Source>(ap.get(), timePlace.get());
    double *tmpGrid = ap->getMasterGridJ2000();
    qDebug() << "Source Ra:" << tmpGrid[0] * 180.0/M_PI << "Dec:"
             << tmpGrid[1] * 180.0/M_PI;
    // make a telescope
    telescope = std::make_unique<Telescope>(ap.get(), timePlace.get(), source.get());
    // get a pointer to the good detectors in the array
    array->updateDetectorIndices();
    int* di = array->getDetectorIndices();

    //setting up physical coordinate system
    // project telescope pointing to mastergrid tangent
    telescope->absToPhysEqPointing();
    // generate pointing for each detector
    for (int i = 0; i < array->getNDetectors(); ++i) {
        array->detectors[di[i]].getPointing(telescope.get(), timePlace.get(), source.get());
        array->detectors[di[i]].getAzElPointing(telescope.get());
    }
    // find map bounds
    array->findMinMaxXY();

    //despiking at detector level, detector by detector
    // first round of despiking
    for (int i = 0; i < array->getNDetectors(); ++i) {
        array->detectors[di[i]].despike(ap->getDespikeSigma());
    }
    array->updateDetectorIndices();
    di = array->getDetectorIndices();

    //replace flagged data in scans with faked data (useful for pca if needed)
    array->fakeFlaggedData(telescope.get());

    //make a fake source for psf determination
    for(int i = 0; i <array->getNDetectors(); ++i){
        array->detectors[di[i]].makeKernelTimestream(telescope.get());
    }

    //lowpass the data
    for(int i = 0; i < array->getNDetectors(); ++i){
        array->detectors[di[i]].lowpass(&array->digFiltTerms[0], array->nFiltTerms);
    }
    //clean out overflagged scans
    VecBool obsFlags(array->detectors[di[0]].getNSamples());
    for (int j = 0; j < array->detectors[di[0]].getNSamples(); ++j) {
        obsFlags[j] = 0;
        for(int i = 0; i < array->getNDetectors(); ++i)
            if (array->detectors[di[i]].hSampleFlags[j]) obsFlags[j] = 1;
    }
    telescope->checkScanLengths(obsFlags, array->detectors[0].getSamplerate());

    //generate the pointing signals
    // estimate average extinction
    for(int i = 0;i < array->getNDetectors(); ++i){
        array->detectors[di[i]].estimateExtinction(array->getAvgTau());
    }

    //create vectors to store future fit parameters for iterative cleaning
    fitParams.assign(array->getNDetectors(), 7, 0.);
    previousFitParams.assign(array->getNDetectors(), 7, 1.0);

    //create the vector that stores whether to terminate iteration for each detector
    needsIteration.assign(array->getNDetectors(), 1.0);

    //store original hValues for iterative cleaning
    int nSamples = array->detectors[di[0]].getNSamples();
    originalhVals.assign(array->getNDetectors(), nSamples, 0.);
    for(int i = 0; i < array->getNDetectors(); ++i){
        for(int j = 0; j < nSamples; j++){
            originalhVals[i][j] = array->detectors[di[i]].hValues[j];
        }
    }

    //grab the iteration loop cap and the percent change cutoff
    iteration = 0;
    cap = ap->getCleanIterationCap();
    cutOff = ap->getCleanIterationCutoff();

    // generate original map before cleaning
    observation = std::make_unique<Observation>(ap.get());
    observation->generateBeammaps(array.get(), telescope.get());
    showMap(0);
}

void MainWindow::showMap(int did)
{
    int ncols = observation->ncols;
    int nrows = observation->nrows;
    colorMap->data()->setSize(ncols, nrows);
    colorMap->data()->setRange(QCPRange(0, 2), QCPRange(0, 2));
    double z;
    double vmin = 1e10, vmax = -1e10;
    for (int x=0; x < ncols; ++x)
    {
        for (int y=0; y < nrows; ++y)
        {
            z = observation->beammapSignal[did][y * ncols + x];
            vmin = (z != 0. && z < vmin)?z:vmin;
            vmax = (z != 0. && z > vmax)?z:vmax;
            colorMap->data()->setCell(x, y, (z==0.)?NAN:z);
        }
    }
    colorMap->setGradient(QCPColorGradient::gpPolar);
    colorMap->setDataRange(QCPRange(vmin, vmax));
    // colorMap->rescaleDataRange(true);
    ui->canvas->rescaleAxes();
    ui->canvas->replot();
}

void MainWindow::showStatus(QMouseEvent* event)
{
    double x = ui->canvas->xAxis->pixelToCoord(event->pos().x());
    double y = ui->canvas->yAxis->pixelToCoord(event->pos().y());
    int i = 0;
    int j = 0;
    colorMap->data()->coordToCell(x, y, &j, &i);
    double z = colorMap->data()->cell(j, i);
    QString tooltip = QString("x=%1 y=%2 z=%3").arg(x).arg(y).arg(z);
    ui->statusBar->showMessage(tooltip);
}

void MainWindow::runBeammap()
{
    /*
    int nFiles = ap->getNFiles();
    qDebug() << "number of ap files: " << nFiles;
    //begin loop over input files
    for(int fileNum = 0; fileNum < nFiles; ++fileNum){
        qDebug() << "starting file number " << fileNum;
        ap->setDataFile(fileNum);
        TimePlace* timePlace = new TimePlace(ap.get());
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
    */
}

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    // read the ap.xml if supplied

    QStringList args = QCoreApplication::arguments();
    if (args.count() > 2)
        qDebug() << "ignore extra commandline arguments:" << args.mid(2);
    if (args.count() > 1)
        w.openApXml(args.at(1));
    return a.exec();
}