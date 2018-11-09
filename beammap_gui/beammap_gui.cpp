#include <QApplication>
#include <QtGlobal>
#include <QDebug>
#include <QShowEvent>
#include <QFileDialog>
#include <QSortFilterProxyModel>
#include <QMessageBox>
#include <QDomDocument>

#include "beammap_gui.h"
#include "dommodel.h"
#include "domitem.h"
#include "ui_beammap_gui.h"
#include "imageview.h"

#include "AnalParams.h"
#include "Observation.h"
#include "TimePlace.h"
#include "Clean.h"
#include "CleanSelector.h"
#include "astron_utilities.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->detComboBox->setModel(&detList);
    qDebug() << "mainwindow initialized";
}

MainWindow::~MainWindow()
{
    delete ui;
    qDebug() << "mainwindow destroyed";
}

void MainWindow::showEvent(QShowEvent* event)
{
    QMainWindow::showEvent(event);
    QApplication::processEvents();
    if(!startUp || event->spontaneous()) return;
    // run statup setting
    startUp = false;
    // read the ap.xml if supplied
    QStringList args = QCoreApplication::arguments();
    if (args.count() > 2)
        qDebug() << "ignore extra commandline arguments:" << args.mid(2);
    if (args.count() > 1)
    {
        QString startUpApXmlPath = args.at(1);
        QMetaObject::invokeMethod(this, [this, startUpApXmlPath] {
            this->openApXml(startUpApXmlPath);
        }, Qt::ConnectionType::QueuedConnection);
    }
    return;
}

void MainWindow::openApXml()
{
     QString path = QFileDialog::getOpenFileName(
                this, tr("Open File"),
                apXmlPath, tr("XML files (*.xml)"));
     QApplication::processEvents();
     if (!path.isEmpty()) {
         QMetaObject::invokeMethod(this, [this, path] {
             openApXml(path);
         }, Qt::ConnectionType::QueuedConnection);
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
            delete apXml; // get rid of the old ap dommodel
            apXml = m;  // and store the new ap dommodel
            apXmlPath = path;  // and store the current path
            // try setup macana
            try {
                qDebug() << "initialize macana core from" << path;
                mAnalParams = std::make_unique<AnalParams>(path.toStdString(), 1);
                qDebug() << "number of observations:" << mAnalParams->getNFiles();
                qDebug() << "macana core initialized";
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
            // set the current observation to the first in the list
            // this will trigger the observation initialization
            ui->obsComboBox->setCurrentIndex(obsId);  // default to the first observation
        }
        fo.close();
    }
}

void MainWindow::setObsId(int idx)
{
    this->obsId = idx;
    // update the observation tree view
    // figure out the model index
    auto pm = static_cast<QSortFilterProxyModel*>(ui->obsComboBox->model());
    QModelIndex parent = pm->match(
        pm->index(0, 0),
        Qt::DisplayRole,
        QVariant::fromValue(QString("observations")),
        1,
        Qt::MatchExactly | Qt::MatchRecursive
        ).at(0);
    QModelIndex mi = pm->mapToSource(
                    pm->index(obsId, 0, parent));
    // set the obs tree view to current obs
    ui->obsTreeView->setRootIndex(mi);
    DomItem* item = static_cast<DomItem*>(mi.internalPointer());
    QString obsName = item->node().firstChildElement().text();
    qDebug() << "select observation file #" << obsId
             << ui->obsComboBox->currentText()
             << obsName;
    // reset detComboBox
    ui->detComboBox->setEnabled(false);
    QApplication::processEvents(); // update the views

    // check macana core
    if (!mAnalParams) {
        qDebug() << "observation is not initialized because macana core is not initialized";
        return;
    }
    // display a message box to show status
    QMessageBox* msg = new QMessageBox(
                QMessageBox::Information,
                "",
                QString("Initializing observation #%1\n\n\t%2").arg(obsId).arg(obsName),
                nullptr, this, Qt::FramelessWindowHint);
    // msg->setStandardButtons(nullptr);  // need to remove the ok button
    connect(this, &MainWindow::observationInitialized, msg, &QMessageBox::accept);
    QMetaObject::invokeMethod(msg, [this] {
        this->initializeObservation();  // initialize macana core objects
    }, Qt::ConnectionType::QueuedConnection);
    if (msg->exec() == QDialog::Accepted) {
        qDebug() << "observation #" << obsId << obsName << "initialized";
        // update the detecor combobox
        updateDetectorIndices();
    } else {
        qDebug() << "failed to initialize observation #" << obsId << obsName << "initialized";
    }
}

void MainWindow::initializeObservation()
{
    qDebug() << "initialize observation #" << obsId;
    mAnalParams->setDataFile(obsId);
    // create an Array object
    mArray = std::make_unique<Array>(mAnalParams.get());
    mArray->populate();
    // make a TimePlace
    mTimePlace = std::make_unique<TimePlace>(mAnalParams.get());
    // make a Source
    // this will set the master grid from the input file
    mSource = std::make_unique<Source>(mAnalParams.get(), mTimePlace.get());
    double *tmpGrid = mAnalParams->getMasterGridJ2000();  // radians
    qDebug() << "Source Ra:" << tmpGrid[0] * DEG_RAD << "Dec:"
             << tmpGrid[1] * DEG_RAD;
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

    //despiking at detector level, detector by detector
    // first round of despiking
    for (int i = 0; i < mArray->getNDetectors(); ++i) {
        mArray->detectors[di[i]].despike(mAnalParams->getDespikeSigma());
    }
    mArray->updateDetectorIndices();
    di = mArray->getDetectorIndices();

    //replace flagged data in scans with faked data (useful for pca if needed)
    mArray->fakeFlaggedData(mTelescope.get());

    //make a fake source for psf determination
    for(int i = 0; i < mArray->getNDetectors(); ++i) {
        mArray->detectors[di[i]].makeKernelTimestream(mTelescope.get());
    }

    //lowpass the data
    for(int i = 0; i < mArray->getNDetectors(); ++i) {
        mArray->detectors[di[i]].lowpass(&mArray->digFiltTerms[0], mArray->nFiltTerms);
    }
    //clean out overflagged scans
    VecBool obsFlags(mArray->detectors[di[0]].getNSamples());
    for (int j = 0; j < mArray->detectors[di[0]].getNSamples(); ++j) {
        obsFlags[j] = 0;
        for(int i = 0; i < mArray->getNDetectors(); ++i)
            if (mArray->detectors[di[i]].hSampleFlags[j]) obsFlags[j] = 1;
    }
    mTelescope->checkScanLengths(obsFlags, mArray->detectors[0].getSamplerate());

    //generate the pointing signals
    // estimate average extinction
    for(int i = 0;i < mArray->getNDetectors(); ++i){
        mArray->detectors[di[i]].estimateExtinction(mArray->getAvgTau());
    }

    //store original hValues for iterative cleaning
    int nSamples = mArray->detectors[di[0]].getNSamples();
    originalhVals.assign(mArray->getNDetectors(), nSamples, 0.);

    for(int i = 0; i < mArray->getNDetectors(); ++i){
        for(int j = 0; j < nSamples; j++){
            originalhVals[i][j] = mArray->detectors[di[i]].hValues[j];
        }
    }
    // end of prep
    emit observationInitialized();
}

/*
    //create vectors to store future fit parameters for iterative cleaning
    fitParams.assign(mArray->getNDetectors(), 7, 0.);
    previousFitParams.assign(mArray->getNDetectors(), 7, 1.0);

    //create the vector that stores whether to terminate iteration for each detector
    needsIteration.assign(mArray->getNDetectors(), 1.0);

    //grab the iteration loop cap and the percent change cutoff
    iteration = 0;
    cap = ap->getCleanIterationCap();
    cutOff = ap->getCleanIterationCutoff();

    // generate original map before cleaning
    observation = std::make_unique<Observation>(ap.get());
    observation->generateBeammaps(array.get(), telescope.get());
    showMap(0);
*/

void MainWindow::setDetId(int idx)
{
    detId = idx;
    foreach (QMdiSubWindow *window, ui->mdiArea->subWindowList()) {
        ImageView *iv = qobject_cast<ImageView*>(window->widget());
        iv->showBeammap(detId);
    }
}

void MainWindow::runBeammap()
{
    qDebug() << "creating beammap of det #" << detId;
    //cleaning
    Clean* cleaner = CleanSelector::getCleaner(mArray.get(), mTelescope.get());
    cleaner->clean();
    delete cleaner;
    updateDetectorIndices();
    if(mArray->getAvgTau() < 0){
        qDebug() << "error: opacity is less than 0.0 on file";
        return;
    }
    // generate beammap Map in imageview after cleanning
    ImageView* iv = new ImageView;
    iv->initializeBeammaps(mAnalParams.get(), mArray.get(), mTelescope.get());
    // show the mdi window
    ui->mdiArea->addSubWindow(iv);
    iv->showMaximized();
    // iv->show();
    // display beammap
    iv->showBeammap(detId);
}

void MainWindow::updateDetectorIndices()
{
    mArray->updateDetectorIndices();
    ui->detComboBox->blockSignals(true); // silence until we finish setup
    QStringList dl;
    int* di = mArray->getDetectorIndices();
    for (int i = 0; i < mArray->getNDetectors(); ++i) {
        dl.push_back(QString("[%1] %2: %3")
                     .arg(i).arg(di[i])
                     .arg(QString::fromStdString(mArray->detectors[di[i]].getName())));
    }
    detList.setStringList(dl);
    ui->detComboBox->setCurrentIndex(detId);
    ui->detComboBox->blockSignals(false);
    ui->detComboBox->setEnabled(true);
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

    /* / load the ap.xml on startup
    auto const c = new QMetaObject::Connection;
    *c = QObject::connect(w, &QMainWindow::, [this, text, c](){
        QObject::disconnect(*c);
        delete c;
    });
    */
    // run the app
    return a.exec();
}