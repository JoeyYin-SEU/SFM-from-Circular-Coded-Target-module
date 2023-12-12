#include "PCL_setting_module.h"

PCL_setting_module::PCL_setting_module(QWidget *parent)
	: QWidget(parent)
{
	ui.setupUi(this);

	connect(ui.lookuptable_color_comboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(setting_cahnge()));
	connect(ui.min_x_doubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(setting_cahnge()));
	connect(ui.max_x_doubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(setting_cahnge()));
	connect(ui.min_y_doubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(setting_cahnge()));
	connect(ui.max_y_doubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(setting_cahnge()));
	connect(ui.min_z_doubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(setting_cahnge()));
	connect(ui.max_z_doubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(setting_cahnge()));
	connect(ui.auto_x_checkbox, SIGNAL(stateChanged(int)), this, SLOT(setting_cahnge()));
	connect(ui.auto_y_checkbox, SIGNAL(stateChanged(int)), this, SLOT(setting_cahnge()));
	connect(ui.auto_z_checkbox, SIGNAL(stateChanged(int)), this, SLOT(setting_cahnge()));
	connect(ui.colorbar_on_checkbox, SIGNAL(stateChanged(int)), this, SLOT(setting_cahnge()));
	connect(ui.grid_on_checkbox, SIGNAL(stateChanged(int)), this, SLOT(setting_cahnge()));
	connect(ui.ori_checkBox, SIGNAL(stateChanged(int)), this, SLOT(setting_cahnge()));
	connect(ui.code_checkBox, SIGNAL(stateChanged(int)), this, SLOT(setting_cahnge()));
	connect(ui.path_checkBox, SIGNAL(stateChanged(int)), this, SLOT(setting_cahnge()));
	connect(ui.grid_on_checkbox, SIGNAL(stateChanged(int)), this, SLOT(setting_cahnge()));
	connect(ui.color_xyz_comboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(setting_cahnge()));
	connect(ui.auto_lut_xyz_checkbox, SIGNAL(stateChanged(int)), this, SLOT(setting_cahnge()));
	connect(ui.min_xyz_lut_doubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(setting_cahnge()));
	
	connect(ui.max_xyz_lut_doubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(setting_cahnge()));
	connect(ui.point_size_spin, SIGNAL(valueChanged(double)), this, SLOT(setting_cahnge()));
	connect(ui.text_scale_spin, SIGNAL(valueChanged(double)), this, SLOT(setting_cahnge()));
	connect(ui.ori_scale_spin, SIGNAL(valueChanged(double)), this, SLOT(setting_cahnge()));
	

	connect(ui.auto_x_checkbox, SIGNAL(stateChanged(int)), this, SLOT(auto_enable_cahnge()));
	connect(ui.auto_y_checkbox, SIGNAL(stateChanged(int)), this, SLOT(auto_enable_cahnge()));
	connect(ui.auto_z_checkbox, SIGNAL(stateChanged(int)), this, SLOT(auto_enable_cahnge()));
	connect(ui.auto_lut_xyz_checkbox, SIGNAL(stateChanged(int)), this, SLOT(auto_enable_cahnge()));

}

PCL_setting_module::~PCL_setting_module()
{}



void PCL_setting_module::setting_cahnge()
{
	if (!in_setting)
	{
		emit setting_cahnge_signal();
	}
}


void PCL_setting_module::auto_enable_cahnge()
{
	ui.max_x_doubleSpinBox->setEnabled(!ui.auto_x_checkbox->isChecked());
	ui.min_x_doubleSpinBox->setEnabled(!ui.auto_x_checkbox->isChecked());
	ui.max_y_doubleSpinBox->setEnabled(!ui.auto_y_checkbox->isChecked());
	ui.min_y_doubleSpinBox->setEnabled(!ui.auto_y_checkbox->isChecked());
	ui.max_z_doubleSpinBox->setEnabled(!ui.auto_z_checkbox->isChecked());
	ui.min_z_doubleSpinBox->setEnabled(!ui.auto_z_checkbox->isChecked());
	ui.max_xyz_lut_doubleSpinBox->setEnabled(!ui.auto_lut_xyz_checkbox->isChecked());
	ui.min_xyz_lut_doubleSpinBox->setEnabled(!ui.auto_lut_xyz_checkbox->isChecked());
}