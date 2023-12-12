#pragma once

#include <QPainter>
#include <QWheelEvent>
#include <iostream>
#include <string>
#include <mutex>
#include <exception>
#include <QLabel>
#include<QPixmap>
#include <QResizeEvent>
class CamShowLabel : public QLabel
{
    Q_OBJECT
public:
    bool can_reset=true;
    CamShowLabel(QWidget* parent);
    ~CamShowLabel();
    void set_image(QImage set_Image);
    void clear_image();
    void OnPresetImage_from_software();
    double cam_ratio;

protected:
    void paintEvent(QPaintEvent* event);
    void mouseDoubleClickEvent(QMouseEvent* event);
    void wheelEvent(QWheelEvent* event) override;
    void mousePressEvent(QMouseEvent* event) override;
    void mouseMoveEvent(QMouseEvent* event) override;
    void mouseReleaseEvent(QMouseEvent* event) override; 

private slots:
    void OnZoomInImage();    
    void OnZoomOutImage();   
    void OnPresetImage();     
    void resizeEvent(QResizeEvent* event);

private:
    std::mutex mtx;
    QImage Image;
    double ZoomValue = 1.0;
    double XPtInterval = 0;
    double YPtInterval = 0;
    QPoint OldPos; 
    bool Pressed = false;
    QString LocalFileName;
};