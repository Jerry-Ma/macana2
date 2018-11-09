#include "imageviewcontrol.h"
#include "ui_imageviewcontrol.h"

ImageViewControl::ImageViewControl(QWidget *parent) :
    QGroupBox(parent),
    ui(new Ui::ImageViewControl)
{
    ui->setupUi(this);
}

ImageViewControl::~ImageViewControl()
{
    delete ui;
}
