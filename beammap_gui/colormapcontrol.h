#ifndef IMAGEVIEWCONTROL_H
#define IMAGEVIEWCONTROL_H

#include <QGroupBox>

namespace Ui {
class ImageViewControl;
}

class ImageViewControl : public QGroupBox
{
    Q_OBJECT

public:
    explicit ImageViewControl(QWidget *parent = nullptr);
    ~ImageViewControl();

private:
    Ui::ImageViewControl *ui;
};

#endif // IMAGEVIEWCONTROL_H
