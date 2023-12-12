#pragma once

#include <QWidget>
#include "ui_PCL_setting_module.h"

class PCL_setting_module : public QWidget
{
	Q_OBJECT

public:
	PCL_setting_module(QWidget *parent = nullptr);
	~PCL_setting_module();
	Ui::PCL_setting_moduleClass ui;
	bool in_setting = false;

private:

private slots:
	void setting_cahnge();
	void auto_enable_cahnge();

signals:
	void setting_cahnge_signal();
};
