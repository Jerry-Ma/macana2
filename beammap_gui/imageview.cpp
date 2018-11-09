#include <QVBoxLayout>
#include "imageview.h"
#include "ui_imageview.h"

#include "Observation.h"

ImageView::ImageView(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::ImageView)
{
    ui->setupUi(this);
    // create maps with a 2x3 grid layout
    // also create the controls alongside
    // create controls
    auto controlLayout = new QVBoxLayout;
    ui->dockWidgetContents->setLayout(controlLayout);
    auto canvasLayout = ui->canvas->plotLayout();
    canvasLayout->clear(); // let's start from scratch and remove the default axis rect
    QList<QCPAxis*> axes;
    for (int n = 0; n < colorMapNames.size(); ++n)
    {
        auto name = colorMapNames[n];
        auto control = new ColorMapControl;
        auto title = name;
        title.replace(0, 1, title[0].toUpper());
        control->setTitle(title);
        controlLayout->addWidget(control);

        QCPAxisRect* ax = new QCPAxisRect(ui->canvas);
        axes << ax->axes();
        auto [i, j] = std::div(n, 3);
        canvasLayout->addElement(i, j, ax);

        auto colorMap = new QCPColorMap(ax->axis(QCPAxis::atBottom), ax->axis(QCPAxis::atLeft));
        colorMap->setInterpolate(false);
        colorMap->setGradient(QCPColorGradient::gpPolar);
        control->addColorMap(colorMap);

        colorMapControls[name] = control;
        colorMaps[name] = colorMap;
    }
    // since we've created the axis rects and axes from scratch, we need to place them on
    // according layers, if we don't want the grid to be drawn above the axes etc.
    // place the axis on "axes" layer and grids on the "grid" layer, which is below "axes":
    for (auto ax: axes)
    {
        ax->setLayer("axes");
        ax->grid()->setLayer("grid");
    }
    // QCPColorScale *cbAx0 = new QCPColorScale(ui->canvas);
    // canvasLayout->addElement(0, 3, cbAx0);
    // residualMap->setColorScale(cbAx1);
    ui->canvas->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    // ui->canvas->setSelectionRectMode(QCP::srmZoom);
    // ui->canvas->setOpenGl(true);
    // lock the signal model and residual control
    colorMapControls["model"]->hide();
    colorMapControls["residual"]->hide();
}

ImageView::~ImageView()
{
    delete ui;
}

void ImageView::initializeBeammaps(AnalParams* ap, Array* ar, Telescope* tel)
{
    mAnalParams = ap;
    mArray = ar;
    mTelescope = tel;
    mObservation = std::make_unique<Observation>(mAnalParams);
    mObservation->generateBeammaps(mArray, mTelescope);
}

void ImageView::showBeammap(int detId)
{
    qDebug() << "plot beammap of detector #" << detId;

    auto fixZero = [] (double z) -> double {
        return (z == 0.)? std::nan(""): z;
    };
    setColorMapData("signal", mObservation->beammapSignal[detId], fixZero);
    setColorMapData("sigma", mObservation->beammapSigma[detId], fixZero);
    setColorMapData("weight", mObservation->beammapWeight[detId], fixZero);
}

void ImageView::setColorMapData(QString name, double* data, std::function<double(double)> func)
{
    int ncols = mObservation->ncols;
    int nrows = mObservation->nrows;
    colorMaps[name]->data()->setSize(ncols, nrows);
    colorMaps[name]->data()->setRange(QCPRange(0, 2), QCPRange(0, 2));
    double z;
    for (int x=0; x < ncols; ++x)
        for (int y=0; y < nrows; ++y)
        {
            z = func(data[y * ncols + x]);
            colorMaps[name]->data()->setCell(x, y, z);
        }
    // remove the cache so that it can be rebuilt later
    ImageView::colorMapDataCache.remove(colorMaps[name]);
    // stat the data to by displayed in quantile. This will rebuild the cache
    colorMapControls[name]->setCutLevelQuantile(0.005);
}

void ImageView::showStatus(QMouseEvent* event)
{
    // find out with axis the mouse is at
    int xpos = event->pos().x();
    int ypos = event->pos().y();
    QString currentName;
    for (auto name: colorMapNames)
    {
        if (colorMaps[name]->keyAxis()->axisRect()->rect().contains(event->pos()))
        {
            currentName = name;
            break;
        }
    }
    if (currentName.isEmpty()) return;
    double x = colorMaps[currentName]->keyAxis()->pixelToCoord(xpos);
    double y = colorMaps[currentName]->valueAxis()->pixelToCoord(ypos);
    int i = 0;
    int j = 0;
    colorMaps[currentName]->data()->coordToCell(x, y, &j, &i);
    double z = colorMaps[currentName]->data()->cell(j, i);
    QString tooltip = QString("%1: x=%2 y=%3 z=%4").arg(currentName).arg(x).arg(y).arg(z);
    ui->statusBar->showMessage(tooltip);
}

const std::vector<double>& ImageView::getOrCreateColorMapDataCache(QCPColorMap* map)
{
    if (!ImageView::colorMapDataCache.contains(map))
    {
        int ncols = map->data()->keySize();
        int nrows = map->data()->valueSize();
        double z;
        ImageView::colorMapDataCache[map].clear();
        ImageView::colorMapDataCache[map].reserve(static_cast<size_t>(ncols * nrows));
        for (int x = 0; x < ncols; ++x)
            for (int y = 0; y < nrows; ++y)
            {
                z = map->data()->cell(x, y);
                if (!std::isnan(z) && z != 0.)
                    ImageView::colorMapDataCache[map].push_back(z);
            }
        std::sort(ImageView::colorMapDataCache[map].begin(), ImageView::colorMapDataCache[map].end());
        qDebug() << "cached" << ImageView::colorMapDataCache[map].size() << "of" << ncols * nrows << "good data points";
    }
    return ImageView::colorMapDataCache[map];
}

QMap<QCPColorMap*, std::vector<double>> ImageView::colorMapDataCache;