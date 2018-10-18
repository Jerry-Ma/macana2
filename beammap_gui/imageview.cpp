#include "imageview.h"
#include "ui_imageview.h"

ImageView::ImageView(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::ImageView)
{
    ui->setupUi(this);
}

ImageView::~ImageView()
{
    delete ui;
}
