#ifndef BEAMMAP_GUI_H
#define BEAMMAP_GUI_H

#include <QMainWindow>

class DomModel;

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
    void runBeammap();

private:
    Ui::MainWindow *ui;
    DomModel *model;
    QString xmlPath;
};

#endif // BEAMMAP_GUI_H
