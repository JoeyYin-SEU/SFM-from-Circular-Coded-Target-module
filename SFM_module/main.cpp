#include "SFM_module.h"
#include <QtWidgets/QApplication>
#include <stdlib.h>
#include <Windows.h>
#include <QStyleFactory>

int main(int argc, char* argv[])
{

    HWND hwnd;
    hwnd = FindWindow(L"ConsoleWindowClass", NULL);
    if (hwnd)
    {
        ShowWindow(hwnd, SW_HIDE);
    }

    QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
    QCoreApplication::setAttribute(Qt::AA_UseHighDpiPixmaps);
    QGuiApplication::setHighDpiScaleFactorRoundingPolicy(Qt::HighDpiScaleFactorRoundingPolicy::PassThrough);
    QApplication a(argc, argv);
    a.setOrganizationName("@ZhuoyiYin");
    a.setApplicationName("Photomechanics_Engine_Tools_CameraCalibration");
    a.setStyle(QStyleFactory::create("fusion"));
    SFM_module w;
    w.show();
    return a.exec();
}
