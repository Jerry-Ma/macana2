#ifndef IMAGEVIEWCONTROL_H
#define IMAGEVIEWCONTROL_H

#include <QGroupBox>
#include <QMap>
#include "qcustomplot/qcustomplot.h"



#include <QAbstractItemModel>
#include <QModelIndex>
#include <QVariant>

/* This class stores settings for the color map
 * It provides a model based interface to be manipulated by multiple view widgets
 */
class ColorMapSettings : public QAbstractItemModel
{
    Q_OBJECT

public:
    ColorMapSettings(QObject* parent = nullptr);
    ~ColorMapSettings() override;

    QVariant data(const QModelIndex &index, int role) const override;
    QVariant headerData(int section, Qt::Orientation orientation,
                        int role = Qt::DisplayRole) const override;

    QModelIndex index(int row, int column,
                      const QModelIndex &parent = QModelIndex()) const override;
    QModelIndex parent(const QModelIndex &index) const override;

    int rowCount(const QModelIndex &parent = QModelIndex()) const override;
    int columnCount(const QModelIndex &parent = QModelIndex()) const override;

    Qt::ItemFlags flags(const QModelIndex &index) const override;
    bool setData(const QModelIndex &index, const QVariant &value,
                 int role = Qt::EditRole) override;

    // setter going through model
    void set(const QString& key, const QVariant& value);
    // template getter
    template<typename T>
    const T get(const QString& key) const
    {
        return (*this)[key].value<T>();
    }

    // low level getters and setters
    // directly to the underlaying data
    QVariant& operator[](const QString& key);
    const QVariant operator[](const QString& key) const;
    QVariant& operator[](int i);
    const QVariant operator[](int i) const;

    // key tricks
    const QString keyOf(int i) const;
    int indexOf(const QString& key);

signals:
    void cutLevelChanged(double, double);
    void rangeChanged(double, double);

public slots:
    void setCutLevel(double, double);
    void setRange(double, double);

private:
    QMap<QString, QVariant> mSettings {
        {"name", ""},
        {"minValue", 0.},
        {"maxValue", 1.},
        {"lowerCutLevel", 0.},
        {"upperCutLevel", 1.},
    };
};


namespace Ui {
class ColorMapControl;
}

class ColorMapControl : public QGroupBox
{
    Q_OBJECT

public:
    explicit ColorMapControl(QWidget *parent = nullptr);
    ~ColorMapControl();

    void addColorMap(QCPColorMap* colorMap, QString name="", bool asMaster=false);
    void removeColorMap(QCPColorMap* colorMap);
    void clearColorMap();
    void setMasterColorMap(QCPColorMap* colorMap);
    void setCutLevelQuantile(double quantile=0.005);

private:
    Ui::ColorMapControl *ui;

    // maps to be managed
    QMap<QCPColorMap*, QString> colorMaps;
    // the master is responsible for the settings of all the maps
    QCPColorMap* masterColorMap;

    ColorMapSettings settings;
    QDataWidgetMapper* settingsMapper;
    void applyCutLevels();
};

#endif // IMAGEVIEWCONTROL_H
