#include "SFM_setting_module.h"

SFM_setting_module::SFM_setting_module(QWidget *parent)
	: QWidget(parent)
{
	ui.setupUi(this);

	connect(ui.r_k_doubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(setting_cahnge()));
	connect(ui.r_k_1_doubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(setting_cahnge()));
	connect(ui.r_k_2_doubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(setting_cahnge()));
	connect(ui.min_radius_doubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(setting_cahnge()));
	connect(ui.max_radius_doubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(setting_cahnge()));
	connect(ui.max_aspect_ratio_doubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(setting_cahnge()));
	connect(ui.max_radius_error_doubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(setting_cahnge()));
	connect(ui.min_contour_number_spinBox, SIGNAL(valueChanged(int)), this, SLOT(setting_cahnge()));
	connect(ui.min_contour_points_spinBox, SIGNAL(valueChanged(int)), this, SLOT(setting_cahnge()));
	connect(ui.max_foreground_std_spinBox, SIGNAL(valueChanged(int)), this, SLOT(setting_cahnge()));
	connect(ui.max_background_std_spinBox, SIGNAL(valueChanged(int)), this, SLOT(setting_cahnge()));
	
	
	connect(ui.min_gap_spinBox, SIGNAL(valueChanged(int)), this, SLOT(setting_cahnge()));
	//connect(ui.min_view_spinBox, SIGNAL(valueChanged(int)), this, SLOT(setting_cahnge()));
	connect(ui.color_comboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(setting_cahnge()));
	connect(ui.code_bits_comboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(setting_cahnge()));
	connect(ui.feature_detect_method_comboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(setting_cahnge()));
	connect(ui.sub_pixel_method_comboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(setting_cahnge()));
	

	connect(ui.optimize_checkBox, SIGNAL(stateChanged(int)), this, SLOT(using_fixed_cxcy()));
	connect(ui.fixed_camera_checkBox, SIGNAL(stateChanged(int)), this, SLOT(using_fixed_all()));

}

SFM_setting_module::~SFM_setting_module()
{}


void SFM_setting_module::setting_cahnge()
{
	emit setting_cahnge_signal();
}

void SFM_setting_module::using_fixed_cxcy()
{
	if (ui.optimize_checkBox->isChecked())
	{
		ui.min_cxcy_BA_spinBox->setEnabled(false);
	}
	else
	{
		ui.min_cxcy_BA_spinBox->setEnabled(true);
	}
}

void SFM_setting_module::set_range()
{

	if (ui.unlimit_checkBox->isChecked())
	{
		ui.limit_spinBox->setEnabled(false);
	}
	else
	{
		ui.limit_spinBox->setEnabled(true);
	}
}

void SFM_setting_module::using_fixed_all()
{
	if (ui.fixed_camera_checkBox->isChecked())
	{
		ui.fixed_f_checkBox->setEnabled(false);
		ui.fixed_c_checkBox->setEnabled(false);
		ui.optimize_checkBox->setEnabled(false);
		ui.min_cxcy_BA_spinBox->setEnabled(false);
		ui.skew_checkBox->setEnabled(false);
		ui.uniform_f_checkBox->setEnabled(false);
		ui.unlimit_checkBox->setEnabled(false);
		ui.k_distort_comboBox->setEnabled(false);
		ui.p_distort_comboBox->setEnabled(false);
		ui.s_distort_comboBox->setEnabled(false);
		ui.limit_spinBox->setEnabled(false);
	}
	else
	{
		ui.fixed_f_checkBox->setEnabled(true);
		ui.fixed_c_checkBox->setEnabled(true);
		ui.optimize_checkBox->setEnabled(true);
		ui.min_cxcy_BA_spinBox->setEnabled(true);
		ui.skew_checkBox->setEnabled(true);
		ui.uniform_f_checkBox->setEnabled(true);
		ui.unlimit_checkBox->setEnabled(true);
		ui.k_distort_comboBox->setEnabled(true);
		ui.p_distort_comboBox->setEnabled(true);
		ui.s_distort_comboBox->setEnabled(true);
		ui.limit_spinBox->setEnabled(true);
	}
}