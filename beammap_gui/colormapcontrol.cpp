#include "colormapcontrol.h"
#include "ui_colormapcontrol.h"

#include "imageview.h"
#include "mathutils.h"


ColorMapSettings::ColorMapSettings(QObject *parent)
    : QAbstractItemModel(parent)
{
}
ColorMapSettings::~ColorMapSettings() {}

int ColorMapSettings::columnCount(const QModelIndex & /* parent */) const
{
    return mSettings.size();
}

QVariant ColorMapSettings::data(const QModelIndex &index, int role) const
{
    if (!index.isValid())
        return QVariant();

    if (role != Qt::DisplayRole && role != Qt::EditRole)
        return QVariant();

    return (*this)[index.column()];
}

Qt::ItemFlags ColorMapSettings::flags(const QModelIndex &index) const
{
    if (!index.isValid())
        return nullptr;
    return Qt::ItemIsEditable | QAbstractItemModel::flags(index);
}

QVariant ColorMapSettings::headerData(int section, Qt::Orientation orientation,
                               int role) const
{
    if (orientation == Qt::Horizontal && role == Qt::DisplayRole)
    {
        return keyOf(section);
    }
    return QVariant();
}

QModelIndex ColorMapSettings::index(int row, int column, const QModelIndex &parent) const
{
    if (parent.isValid() && parent.column() != 0) return QModelIndex();
    return createIndex(row, column, nullptr);
}

QModelIndex ColorMapSettings::parent(const QModelIndex &index) const
{
    Q_UNUSED(index);
    return QModelIndex();
}

int ColorMapSettings::rowCount(const QModelIndex &parent) const
{
    Q_UNUSED(parent);
    return 1;
}

bool ColorMapSettings::setData(const QModelIndex &index, const QVariant &value, int role)
{
    if (role != Qt::EditRole) return false;
    (*this)[index.column()] = value;
    emit dataChanged(index, index);
    auto key = keyOf(index.column());
    if (key.endsWith("CutLevel"))
    {
        emit cutLevelChanged(get<double>("lowerCutLevel"), get<double>("upperCutLevel"));
    } else if (key == "minValue" || key == "maxValue")
    {
        emit rangeChanged(get<double>("minValue"), get<double>("maxValue"));
    }
    return true;
}

void ColorMapSettings::set(const QString& key, const QVariant& value)
{
    auto i = indexOf(key);
    setData(index(0, i, QModelIndex()), value);
}

void ColorMapSettings::setCutLevel(double lower, double upper)
{
    set("lowerCutLevel", lower);
    set("upperCutLevel", upper);
}


void ColorMapSettings::setRange(double min, double max)
{
    set("minValue", min);
    set("maxValue", max);
}




int ColorMapSettings::indexOf(const QString& key)
{
    return mSettings.keys().indexOf(key);
}

const QString ColorMapSettings::keyOf(int i) const
{
    auto it = mSettings.constBegin();
    it += i;
    return it.key();
}

QVariant& ColorMapSettings::operator[](const QString& key)
{
    return mSettings[key];
}

const QVariant ColorMapSettings::operator[](const QString& key) const
{
    return mSettings[key];
}

QVariant& ColorMapSettings::operator[](int i)
{
    return mSettings[keyOf(i)];
}

const QVariant ColorMapSettings::operator[](int i) const
{
    return mSettings[keyOf(i)];
}

ColorMapControl::ColorMapControl(QWidget *parent) :
    QGroupBox(parent),
    ui(new Ui::ColorMapControl)
{
    ui->setupUi(this);
    settingsMapper = new QDataWidgetMapper(this);
    settingsMapper->setModel(&settings);
    settingsMapper->addMapping(ui->lowerCutLevelSpinBox, settings.indexOf("lowerCutLevel"));
    settingsMapper->addMapping(ui->upperCutLevelSpinBox, settings.indexOf("upperCutLevel"));
    settingsMapper->addMapping(ui->cutLevelSlider, settings.indexOf("lowerCutLevel"), "lowerValue");
    settingsMapper->addMapping(ui->cutLevelSlider, settings.indexOf("upperCutLevel"), "upperValue");
    settingsMapper->addMapping(ui->cutLevelSlider, settings.indexOf("minValue"), "minimumValue");
    settingsMapper->addMapping(ui->cutLevelSlider, settings.indexOf("maxValue"), "maximumValue");
    settingsMapper->toFirst();

    // update the plot
    connect(&settings, &ColorMapSettings::cutLevelChanged, this, &ColorMapControl::applyCutLevels);

    // interact with the double slider
    connect(&settings, &ColorMapSettings::cutLevelChanged, ui->cutLevelSlider, [this] (double lower, double upper) {
        this->ui->cutLevelSlider->blockSignals(true);
        this->ui->cutLevelSlider->setValues(lower, upper);
        this->ui->cutLevelSlider->blockSignals(false);
    });
     connect(&settings, &ColorMapSettings::rangeChanged, ui->cutLevelSlider, [this] (double min, double max) {
        this->ui->cutLevelSlider->blockSignals(true);
        this->ui->cutLevelSlider->setRange(min, max);
        this->ui->cutLevelSlider->blockSignals(false);
    });
    connect(ui->cutLevelSlider, &DoubleRangeSlider::valuesChanged, &settings, &ColorMapSettings::setCutLevel);
    connect(ui->cutLevelSlider, &DoubleRangeSlider::rangeChanged, &settings, &ColorMapSettings::setRange);
}

ColorMapControl::~ColorMapControl()
{
    delete ui;
}


void ColorMapControl::addColorMap(QCPColorMap* colorMap, QString name, bool asMaster)
{
    colorMaps[colorMap] = name;
    if (colorMaps.size() == 1) asMaster = true;
    if (asMaster) setMasterColorMap(colorMap);
}

void ColorMapControl::removeColorMap(QCPColorMap* colorMap)
{
    colorMaps.remove(colorMap);
    if (colorMap == masterColorMap)
    {
        if (colorMaps.empty())
        {
            setMasterColorMap(nullptr);
        } else {
            setMasterColorMap(colorMaps.firstKey());
        }
    }
}

void ColorMapControl::clearColorMap()
{
    for (auto colorMap: colorMaps.keys())
    {
        removeColorMap(colorMap);
    }
}

void ColorMapControl::setMasterColorMap(QCPColorMap* colorMap)
{
    masterColorMap = colorMap;
}

void ColorMapControl::setCutLevelQuantile(double quantile)
{
    double symQuantile = 1. - quantile;
    if (quantile > symQuantile) std::swap(quantile, symQuantile);

    auto qs = Quantile(ImageView::getOrCreateColorMapDataCache(masterColorMap), {0., quantile, symQuantile, 1.});
    qDebug() << "qs" << qs;
    settings.set("minValue", qs[0]);
    settings.set("lowerCutLevel", qs[1]);
    settings.set("upperCutLevel", qs[2]);
    settings.set("maxValue", qs[3]);
}

void ColorMapControl::applyCutLevels()
{
    if (colorMaps.empty()) return;
    double lower = settings.get<double>("lowerCutLevel");
    double upper = settings.get<double>("upperCutLevel");
    for (auto colorMap: colorMaps.keys())
    {
        colorMap->setDataRange(QCPRange(lower, upper));
        // colorMap->rescaleDataRange(true);
        colorMap->parentPlot()->rescaleAxes();
        colorMap->parentPlot()->replot();
    }
}
