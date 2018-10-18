#include <QApplication>
#include <QtGlobal>
#include <QDomDocument>

#include "beammap_gui.h"
#include "dommodel.h"
#include "domitem.h"
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
    qDebug() << "mainwindow initialized";
}

MainWindow::~MainWindow()
{
    delete ui;
    qDebug() << "mainwindow destroyed";
}

void MainWindow::openApXml()
{
     QString path = QFileDialog::getOpenFileName(
                this, tr("Open File"),
                apXmlPath, tr("XML files (*.xml)"));
     if (!path.isEmpty()) {
         return openApXml(path);
     }
     return;
}

void MainWindow::openApXml(const QString& path)
{
    qDebug() << "read analysis parameters from" << path;
    QFile fo(path);
    if (fo.open(QIODevice::ReadWrite)) {
        QDomDocument d;
        if (d.setContent(&fo)) {
            auto m = new DomModel(d, this);
            // view entire tree
            ui->apXmlTreeView->setModel(m);
            ui->apXmlTreeView->expandAll();
            // view current observation info
            ui->obsTreeView->setModel(m);
            ui->obsTreeView->setRootIndex(QModelIndex());
            ui->obsTreeView->expandAll();
            // view observation list
            ui->obsComboBox->blockSignals(true); // silence until we finish setup
            // set up a proxy to populate the observation combobox
            QSortFilterProxyModel* p = new QSortFilterProxyModel(this);
            p->setSourceModel(m);
            // find out and keep only the file entries
            p->setRecursiveFilteringEnabled(true);
            p->setFilterRegExp(QRegExp(
                        "f\\d+",
                        Qt::CaseInsensitive));
            ui->obsComboBox->setModel(p);
            // get the index for observation tag
            QModelIndex obsindex = p->match(
                    p->index(0, 0),
                    Qt::DisplayRole,
                    QVariant::fromValue(QString("observations")),
                    1,
                    Qt::MatchExactly | Qt::MatchRecursive
                    ).at(0);
            ui->obsComboBox->setRootModelIndex(obsindex);
            ui->obsComboBox->blockSignals(false); // let go the signal
            // apply changes to class members
            delete apXml; // get rid of the old ap dommodel
            apXml = m;  // and store the new ap dommodel
            apXmlPath = path;  // and store the current path
            // try setup macana
            try {
                qDebug() << "initialize macana core from" << path;
                mAnalParams = std::make_unique<AnalParams>(path.toStdString(), 1);
                qDebug() << "number of observations:" << mAnalParams->getNFiles();
                qDebug() << "macana core initialized";
                // set the current observation to the first in the list
                // this will trigger the observation initialization
                ui->obsComboBox->setCurrentIndex(obsId);  // default to the first observation
                // restore the views' look
                ui->apXmlTreeView->setPalette(QApplication::palette(ui->apXmlTreeView));
            } catch (const AnalParamsError& error) {
                qDebug() << "failed to initialize macana core:" << error.what();
                // reset the core to null
                mAnalParams.reset(nullptr);
                // taint the view with red
                // QPalette p(ui->apXmlTreeView->palette());
                // p.setColor( QPalette::Base, Qt::red);
                // ui->apXmlTreeView->setPalette(p);
            }
        }
        fo.close();
    }
}

void MainWindow::setObsId(int idx)
{
    this->obsId = idx;
    // update the observation tree view
    // figure out the model index
    auto p = static_cast<QSortFilterProxyModel*>(ui->obsComboBox->model());
    QModelIndex parent = p->match(
        p->index(0, 0),
        Qt::DisplayRole,
        QVariant::fromValue(QString("observations")),
        1,
        Qt::MatchExactly | Qt::MatchRecursive
        ).at(0);
    QModelIndex i = p->mapToSource(
                    p->index(obsId, 0, parent));
    // set the obs tree view to current obs
    ui->obsTreeView->setRootIndex(i);
    DomItem* item = static_cast<DomItem*>(i.internalPointer());
    qDebug() << "select observation file #" << obsId
             << ui->obsComboBox->currentText()
             << item->node().firstChildElement().text();
    initializeObservation();
}


void MainWindow::initializeObservation()
{
    if (!mAnalParams) {
        qDebug() << "observation is not initialized because macana core is not initialized";
        return;
    }
    qDebug() << "initialize observation #" << obsId;
    mAnalParams->setDataFile(obsId);
    // create an Array object
    mArray = std::make_unique<Array>(mAnalParams.get());
    mArray->populate();
    // make a TimePlace
    mTimePlace = std::make_unique<TimePlace>(mAnalParams.get());
    // make a Source
    mSource = std::make_unique<Source>(mAnalParams.get(), mTimePlace.get());
    double *tmpGrid = mAnalParams->getMasterGridJ2000();
    qDebug() << "Source Ra:" << tmpGrid[0] * 180.0/M_PI << "Dec:"
             << tmpGrid[1] * 180.0/M_PI;
    // make a telescope
    mTelescope = std::make_unique<Telescope>(mAnalParams.get(), mTimePlace.get(), mSource.get());
    // get a pointer to the good detectors in the array
    mArray->updateDetectorIndices();
    int* di = mArray->getDetectorIndices();
    //setting up physical coordinate system
    // project telescope pointing to mastergrid tangent
    mTelescope->absToPhysEqPointing();
    // generate pointing for each detector
    for (int i = 0; i < mArray->getNDetectors(); ++i) {
        mArray->detectors[di[i]].getPointing(mTelescope.get(), mTimePlace.get(), mSource.get());
        mArray->detectors[di[i]].getAzElPointing(mTelescope.get());
    }
    // find map bounds
    mArray->findMinMaxXY();

    /*
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
    */
    emit observationInitialized(true);
}

void MainWindow::showMap(int did)
{
    int ncols = mObservation->ncols;
    int nrows = mObservation->nrows;
    colorMap->data()->setSize(ncols, nrows);
    colorMap->data()->setRange(QCPRange(0, 2), QCPRange(0, 2));
    double z;
    double vmin = 1e10, vmax = -1e10;
    for (int x=0; x < ncols; ++x)
    {
        for (int y=0; y < nrows; ++y)
        {
            z = mObservation->beammapSignal[did][y * ncols + x];
            vmin = (z != 0. && z < vmin)?z:vmin;
            vmax = (z != 0. && z > vmax)?z:vmax;
            colorMap->data()->setCell(x, y, (z==0.)?std::nan(""):z);
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
    // setup logging format
    QString pattern = "[%{if-debug}D%{endif}%{if-info}I%{endif}%{if-warning}W%{endif}%{if-critical}C%{endif}%{if-fatal}F%{endif}] %{message}";
    #ifdef QT_DEBUG
        pattern += " (%{function})";
    #endif
    qSetMessagePattern(pattern);
    qDebug() << "CXX version:" << __cplusplus ;

    // create windows
    QApplication a(argc, argv);
    MainWindow w;
    w.show();

    // read the ap.xml if supplied
    QStringList args = QCoreApplication::arguments();
    if (args.count() > 2)
        qDebug() << "ignore extra commandline arguments:" << args.mid(2);
    if (args.count() > 1)
        w.openApXml(args.at(1));

    // run the app
    return a.exec();
}