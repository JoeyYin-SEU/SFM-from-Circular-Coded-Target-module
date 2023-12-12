#include "camshowlabel.h"
CamShowLabel::CamShowLabel(QWidget* parent) :QLabel(parent)
{
    cam_ratio = 1;
    //InitializeCriticalSection(Image);
}
CamShowLabel::~CamShowLabel()
{

}

void CamShowLabel::clear_image()
{
    mtx.lock();
    Image = QImage();
    mtx.unlock();
}
void CamShowLabel::set_image(QImage set_Image)
{
    mtx.lock();
    Image = set_Image;
    mtx.unlock();
}

//QPainter画图
void CamShowLabel::resizeEvent(QResizeEvent* event)
{
    QSize new_size = event->size();
    if (abs(new_size.width() - new_size.height() * cam_ratio) > 10)
    {
        if (new_size.width() < new_size.height() * cam_ratio) {

            new_size.setHeight(new_size.width() / cam_ratio);
            //this->move(0, (old_size.height() - new_size.height()) / 2);
        }
        else {
            new_size.setWidth(new_size.height() * cam_ratio);
            //this->move((old_size.width() - new_size.width()) / 2, 0);
        }
    }

    this->resize(new_size);
    //if (Image.isNull())
    //{
    //    return;
    //}

    //mtx.lock();
    //auto w = (double)Image.width();
    //auto h = (double)Image.height();
    //mtx.unlock();
    //auto ws = (double)this->width();
    //auto hs = (double)this->height();
    //auto qw = w / h * hs;
    //auto qh = ws / (w / h);
    //if (qw < ws*1.05)
    //{
    //    if(qw>=10)
    //        this->setFixedWidth(qw);
    //    else
    //        this->setFixedWidth(10);
    //    this->setMinimumWidth(10);
    //    this->setMaximumWidth(100 * qw);
    //    //this->setFixedHeight(hs);
    //}
    //else
    //{
    //    //this->setFixedWidth(ws);
    //    if (qh >= 10)
    //        this->setFixedHeight(qh);
    //    else
    //        this->setFixedHeight(10);
    //    this->setMinimumHeight(10);
    //    this->setMaximumHeight(100 * qh);
    //}
}
void CamShowLabel::paintEvent(QPaintEvent* event)
{
    if (Image.isNull())
    {
        QPainter painter(this);
        return;
    }
    QLabel::paintEvent(event);
    QPainter painter(this);
    double width = this->width();  
    double height = this->height();
    if (width / 2 + XPtInterval - width / 2 * ZoomValue > 0)
    {
        XPtInterval = -width/2 + width / 2 * ZoomValue;
    }
    if (width / 2 + XPtInterval + width / 2 * ZoomValue < width)
    {
        XPtInterval = width / 2 - width / 2 * ZoomValue;
    }
    if (height / 2 + YPtInterval - height / 2 * ZoomValue > 0)
    {
        YPtInterval = -height / 2 + height / 2 * ZoomValue;
    }
    if (height / 2 + YPtInterval + height / 2 * ZoomValue < height)
    {
        YPtInterval = height / 2 - height / 2 * ZoomValue;
    }
    painter.translate(width / 2 + XPtInterval, height / 2 + YPtInterval);
    painter.scale(ZoomValue, ZoomValue);
    QRect picRect(-width / 2, -height / 2, width, height);
    mtx.lock();
    painter.drawImage(picRect, Image);
    mtx.unlock();

    QPen pen(Qt::green,0.2);
    painter.setPen(pen);
    painter.drawLine(width / 2 - width / 2, 0 - height / 2, width / 2 - width / 2, height - height / 2);
    //绘制纵向线
    painter.drawLine(0 - width / 2, height / 2 - height / 2, width - width / 2, height / 2 - height / 2);
}
//鼠标滚轮滚动
void CamShowLabel::wheelEvent(QWheelEvent* event)
{
    double zoom_now = ZoomValue;
    int value = event->angleDelta().y();
    if (value > 0)
    {
        OnZoomInImage();
        XPtInterval = ZoomValue / zoom_now * XPtInterval - ((double)event->position().x() - (double)this->width() / 2.0) * (ZoomValue / zoom_now - 1);
        YPtInterval = ZoomValue / zoom_now * YPtInterval - ((double)event->position().y() - (double)this->height() / 2.0) * (ZoomValue / zoom_now - 1);
    }
    else
    {
        OnZoomOutImage();
        XPtInterval = ZoomValue / zoom_now * XPtInterval - ((double)event->position().x() - (double)this->width() / 2.0) * (ZoomValue / zoom_now - 1);
        YPtInterval = ZoomValue / zoom_now * YPtInterval - ((double)event->position().y() - (double)this->height() / 2.0) * (ZoomValue / zoom_now - 1);
    }
    update();
}
//鼠标摁下
void CamShowLabel::mousePressEvent(QMouseEvent* event)
{
    if (event->button() != Qt::RightButton)
    {
        return;
    }
    OldPos = event->pos();
    Pressed = true;
}
//鼠标松开
void CamShowLabel::mouseMoveEvent(QMouseEvent* event)
{
    if (!Pressed)
        return QWidget::mouseMoveEvent(event);

    this->setCursor(Qt::SizeAllCursor);
    QPoint pos = event->pos();
    int xPtInterval = pos.x() - OldPos.x();
    int yPtInterval = pos.y() - OldPos.y();
    XPtInterval += xPtInterval;
    YPtInterval += yPtInterval;
    OldPos = pos; 
    update();
}
//鼠标发射事件
void CamShowLabel::mouseReleaseEvent(QMouseEvent*/*event*/)
{
    Pressed = false;
    setCursor(Qt::ArrowCursor);
}
//图片放大
void CamShowLabel::OnZoomInImage()
{
    ZoomValue = ZoomValue * 1.2;
    update();
}
//图片缩小
void CamShowLabel::OnZoomOutImage()
{
    ZoomValue = ZoomValue / 1.2;
    if (ZoomValue < 1)
    {
        ZoomValue = 1;
        return;
    }
    update();
}
//图片还原
void CamShowLabel::OnPresetImage()
{
    ZoomValue = 1.0;
    XPtInterval = 0;
    YPtInterval = 0;
    update();
}

void CamShowLabel::OnPresetImage_from_software()
{
    ZoomValue = 1.0;
    XPtInterval = 0;
    YPtInterval = 0;
    update();
}
void CamShowLabel::mouseDoubleClickEvent(QMouseEvent* event) 
{
    OnPresetImage();
}
