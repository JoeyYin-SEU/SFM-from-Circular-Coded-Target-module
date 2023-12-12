#pragma once

#include <QWidget>
#include "ui_Calibration_Log_module.h"

class Calibration_Log_module : public QWidget
{
	Q_OBJECT

public:
	Calibration_Log_module(QWidget *parent = nullptr);
	~Calibration_Log_module();
	Ui::Calibration_Log_moduleClass ui;

};
