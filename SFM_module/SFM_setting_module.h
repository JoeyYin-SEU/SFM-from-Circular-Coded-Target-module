#pragma once

#include <QWidget>
#include "ui_SFM_setting_module.h"

class SFM_setting_module : public QWidget
{
	Q_OBJECT

public:
	SFM_setting_module(QWidget *parent = nullptr);
	~SFM_setting_module();
	Ui::SFM_setting_moduleClass ui;

private:

private slots:
	void setting_cahnge();
	void using_fixed_cxcy();
	void using_fixed_all();
	void set_range();
signals:
	void setting_cahnge_signal();
};
