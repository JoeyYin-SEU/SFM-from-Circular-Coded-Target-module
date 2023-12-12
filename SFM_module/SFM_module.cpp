#include "SFM_module.h"

SFM_module::SFM_module(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);
	init_box();
	connect(ui.camera_model_pushButton, SIGNAL(clicked()), this, SLOT(show_setting_window()));
	connect(ui.open_image_pushButton, SIGNAL(clicked()), this, SLOT(open_images()));
	connect(ui.clear_pushButton, SIGNAL(clicked()), this, SLOT(clear_all_data()));
	connect(ui.calculate_pushButton, SIGNAL(clicked()), this, SLOT(calculate_multiview()));
	connect(ui.show_detect_horizontalSlider, SIGNAL(valueChanged(int)), this, SLOT(update_show_view()));
	connect(ui.show_detect_spinBox, SIGNAL(valueChanged(int)), this, SLOT(update_show_view()));
	connect(ui.show_detect_horizontalSlider, SIGNAL(valueChanged(int)), ui.show_detect_spinBox, SLOT(setValue(int)));
	connect(ui.show_detect_spinBox, SIGNAL(valueChanged(int)), ui.show_detect_horizontalSlider, SLOT(setValue(int)));
	connect(ui.actionKeyPoints, SIGNAL(triggered()), this, SLOT(select_path_read_keypoint()));
	connect(ui.actionResults, SIGNAL(triggered()), this, SLOT(select_path_read_result()));
	
	connect(ui.tableWidget, SIGNAL(cellChanged(int, int)), this, SLOT(update_checked_data(int, int)));
	connect(ui.levelling_tableWidget, SIGNAL(cellChanged(int, int)), this, SLOT(update_levelling_data(int, int)));
	connect(ui.export_pushButton, SIGNAL(clicked()), this, SLOT(export_inf_file()));
	connect(ui.save_focal_pixel_pushButton, SIGNAL(clicked()), this, SLOT(save_camera_params()));
	connect(ui.delete_focal_pixel_pushButton, SIGNAL(clicked()), this, SLOT(delete_camera_params()));

	connect(ui.draw_setting_pushButton, SIGNAL(clicked()), this, SLOT(show_draw_setting_window()));
	connect(ui.add_scale_pushButton, SIGNAL(clicked()), this, SLOT(add_scale_lable()));
	connect(ui.clear_scale_pushButton, SIGNAL(clicked()), this, SLOT(clear_scale_lable()));

	connect(ui.reproject_error_limit_checkBox, SIGNAL(stateChanged(int)), this, SLOT(using_reporject_error_limit()));
	connect(ui.reproject_error_limit_doubleSpinBox, SIGNAL(valueChanged(double)), this, SLOT(update_table_view()));

	connect(ui.focal_pixel_comboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(update_camera_params_enable()));
	
}

SFM_module::~SFM_module()
{}
//sbsbsb
std::vector<QString> knuth_sampling(QDir pNum, int m, int n)
{
	std::vector<QString> result;
	for (int i = 0; i < n; i++)
	{
		if (rand() % (n - i) < m)
		{
			result.push_back("C:\\Users\\Yzy_seu\\Desktop\\Image\\" + pNum[i]);
			m--;
		}
	}
	return result;
}
void SFM_module::init_box()
{
	QStringList sl;
	sl << "cpu" << "get" << "NumberOfLogicalProcessors";
	QProcess p;
	p.setProcessChannelMode(QProcess::MergedChannels);
	p.start("wmic", sl);
	p.waitForFinished();
	QString result = QString::fromLocal8Bit(p.readAllStandardOutput());
	result = result.remove(sl.last(), Qt::CaseInsensitive);
	result = result.replace("\r", "");
	result = result.replace("\n", "");
	result = result.simplified();

	Levelling_R = Eigen::Matrix3d::Zero(3, 3);
	Levelling_R(0, 0) = 1;
	Levelling_R(1, 1) = 1;
	Levelling_R(2, 2) = 1;
	Levelling_T = Eigen::Vector3d(0, 0, 0);

	update_camera_params();

	SFM_setting_window = new SFM_setting_module();
	connect(SFM_setting_window, SIGNAL(setting_cahnge_signal()), this, SLOT(setting_window_change()));
	SFM_setting_window->show();
	qApp->processEvents();
	SFM_setting_window->hide();

	cam_log_ui = new Calibration_Log_module();
	cam_log_ui->show();
	cam_log_ui->setWindowOpacity(1);
	qApp->processEvents();
	cam_log_ui->hide();

	PCL_setting_window = new PCL_setting_module();
	connect(PCL_setting_window, SIGNAL(setting_cahnge_signal()), this, SLOT(pcl_setting_window_change()));
	PCL_setting_window->show();
	qApp->processEvents();
	PCL_setting_window->hide();

	bool ok;
	if (result.toInt(&ok) && ok && result.toInt() > 1)
	{
		SFM_setting_window->ui.thread_spinBox->setValue(result.toInt());
	}

	vtkObject::GlobalWarningDisplayOff();

	lut_for_vtk = vtkSmartPointer<vtkLookupTable>::New();
	VTK_plot::set_lut_mode(lut_for_vtk, 0);

	viewer_for_vtk.reset(new pcl::visualization::PCLVisualizer("viewer", false));

	window_for_vtk->AddRenderer(viewer_for_vtk->getRendererCollection()->GetFirstRenderer());
	ui.VTK_xyz_widget->setRenderWindow(window_for_vtk.Get());
	viewer_for_vtk->setupInteractor(ui.VTK_xyz_widget->interactor(), ui.VTK_xyz_widget->renderWindow());
	viewer_for_vtk->setBackgroundColor(1, 1, 1);
	ui.VTK_xyz_widget->update();

	cloud_for_vtk.reset(new pcl::PointCloud<pcl::PointXYZRGBA>);
	cloud_for_vtk->points.resize(50000);
	for (std::size_t ii = 0; ii < cloud_for_vtk->points.size(); ++ii)
	{
		cloud_for_vtk->points[ii].x = 2 * rand() / (RAND_MAX + 1.0f) - 1;
		cloud_for_vtk->points[ii].y = 2 * rand() / (RAND_MAX + 1.0f) - 1;
		cloud_for_vtk->points[ii].z = sqrt(cloud_for_vtk->points[ii].x* cloud_for_vtk->points[ii].x + cloud_for_vtk->points[ii].y * cloud_for_vtk->points[ii].y);
	}
	Points_xyz_range range_now = VTK_plot::get_xyz_range(cloud_for_vtk);

	cubeaxes_for_vtk = vtkSmartPointer<vtkCubeAxesActor>::New();
	VTK_plot::set_grid_mode(viewer_for_vtk, cubeaxes_for_vtk, range_now, true);
	viewer_for_vtk->getRendererCollection()->GetFirstRenderer()->GetActors();
    viewer_for_vtk->getRendererCollection()->GetFirstRenderer()->ResetCamera();
    viewer_for_vtk->getRendererCollection()->GetFirstRenderer()->ResetCameraClippingRange();
    viewer_for_vtk->removeAllPointClouds();
	viewer_for_vtk->removeAllShapes();
	VTK_plot::set_clouds_color_xyz(cloud_for_vtk, lut_for_vtk, range_now, 2);
    viewer_for_vtk->addPointCloud(cloud_for_vtk, "cloud_xyz");
    viewer_for_vtk->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 0.01, "cloud_xyz");

	scalarBarActor_for_vtk = vtkSmartPointer<vtkScalarBarActor>::New();
	scalarbar_Widget_member_for_vtk = vtkSmartPointer<vtkScalarBarWidget>::New();
	scalarBarRep_for_vtk = vtkSmartPointer<vtkScalarBarRepresentation >::New();
	VTK_plot::init_lut(viewer_for_vtk, lut_for_vtk, scalarBarActor_for_vtk, scalarbar_Widget_member_for_vtk, scalarBarRep_for_vtk);
    VTK_plot::set_lut_show(viewer_for_vtk, lut_for_vtk, scalarBarActor_for_vtk, range_now.z_min, range_now.z_max, "Drawing Example");
    viewer_for_vtk->updatePointCloud(cloud_for_vtk, "cloud_xyz");
    viewer_for_vtk->initCameraParameters();
    viewer_for_vtk->resetCamera();
	window_for_vtk->Render();
    viewer_for_vtk->setBackgroundColor(1, 1, 1);
    ui.VTK_xyz_widget->update();
    ui.VTK_xyz_widget->show();


	
}

void SFM_module::using_reporject_error_limit()
{
	ui.reproject_error_limit_doubleSpinBox->setEnabled(ui.reproject_error_limit_checkBox->isChecked());
	update_table_view();
}

void SFM_module::update_checked_data(int rows, int cols)
{
	if (in_reading_keypoint)
	{
		return;
	}
	if (in_setting_keypoint_table)
	{
		return;
	}
	if (rows >= KeyPoint_Enable_serial.size())
	{
		return;
	}
	if (ui.tableWidget->item(rows, cols)->checkState() == Qt::Checked)
	{
		KeyPoint_Enable_serial[rows] = true;
	}
	else if (ui.tableWidget->item(rows, cols)->checkState() == Qt::Unchecked)
	{
		KeyPoint_Enable_serial[rows] = false;
	}
}

void  SFM_module::update_levelling_data(int rows, int cols)
{
	if (in_setting_levelling_table)
	{
		return;
	}
	if (ui.levelling_tableWidget->item(rows, cols)->checkState() == Qt::Checked)
	{
		xyz_need_levelling_serial[rows] = true;
	}
	else if (ui.levelling_tableWidget->item(rows, cols)->checkState() == Qt::Unchecked)
	{
		xyz_need_levelling_serial[rows] = false;
	}
	update_points_view();
}

void SFM_module::update_table_view()
{
	in_setting_keypoint_table = true;
	ui.tableWidget->setRowCount(0);
	ui.tableWidget->setRowCount(Image_serial_name.size());
	std::priority_queue<double> pri_vec_all;
	std::vector<double> vec_all_temp;
	for (int ii = 0; ii < Image_serial_name.size(); ii++)
	{
		if (ui.tableWidget->item(ii, 0) == NULL ||
			(ui.tableWidget->item(ii, 0) &&
				ui.tableWidget->item(ii, 0)->text() == tr("")))
		{
			ui.tableWidget->setItem(ii, 0,
				new QTableWidgetItem(QString("%1").arg(ii + 1, 2, 10, QLatin1Char('0'))));
			if (KeyPoint_Enable_serial[ii])
			{
				ui.tableWidget->item(ii, 0)->setCheckState(Qt::Checked);
			}
			else
			{
				ui.tableWidget->item(ii, 0)->setCheckState(Qt::Unchecked);
			}
			if (ii >= Result_err_every.size()  || !(KeyPoint_Enable_serial[ii]))
			{
				ui.tableWidget->setItem(ii, 1, new QTableWidgetItem(tr("--")));
				ui.tableWidget->setItem(ii, 2, new QTableWidgetItem(tr("--")));
				ui.tableWidget->setItem(ii, 3, new QTableWidgetItem(tr("--")));
				ui.tableWidget->setItem(ii, 4, new QTableWidgetItem(tr("--")));
				QFont nullFont;
				QBrush nullColor(Qt::black);
				ui.tableWidget->item(ii, 1)->setFont(nullFont);
				ui.tableWidget->item(ii, 1)->setForeground(nullColor);
				ui.tableWidget->item(ii, 2)->setFont(nullFont);
				ui.tableWidget->item(ii, 2)->setForeground(nullColor);
				ui.tableWidget->item(ii, 3)->setFont(nullFont);
				ui.tableWidget->item(ii, 3)->setForeground(nullColor);
				ui.tableWidget->item(ii, 4)->setFont(nullFont);
				ui.tableWidget->item(ii, 4)->setForeground(nullColor);
			}
			else
			{
				double mean_err = 0, max_err = -std::numeric_limits<double>::max();
				double median_err = 0;
				int err_number = 0;
				std::priority_queue<double> pri_vec_temp;
				std::vector<double> vec_temp;
				for (int pp = 0; pp < Result_err_every[ii].size(); pp++)
				{
					if (Result_err_every[ii][pp].x == std::numeric_limits<double>::max() || Result_err_every[ii][pp].y == std::numeric_limits<double>::max())
					{
						continue;
					}
					double er_now = sqrt(pow(Result_err_every[ii][pp].x, 2) + pow(Result_err_every[ii][pp].y, 2));
					pri_vec_temp.push(er_now);
					pri_vec_all.push(er_now);
					mean_err += er_now;
					err_number++;
					max_err = max_err < er_now ? er_now : max_err;
				}
				while (pri_vec_temp.size())
				{
					vec_temp.push_back(pri_vec_temp.top());
					pri_vec_temp.pop();
				}
				ui.tableWidget->setItem(ii, 1, new QTableWidgetItem(QString::number(err_number)));
				QFont number_Font;
				QBrush number_Color(Qt::black);
				ui.tableWidget->item(ii, 1)->setFont(number_Font);
				ui.tableWidget->item(ii, 1)->setForeground(number_Color);
				if (err_number == 0)
				{
					ui.tableWidget->item(ii, 0)->setCheckState(Qt::Unchecked);
					ui.tableWidget->setItem(ii, 2, new QTableWidgetItem(tr("NULL")));
					ui.tableWidget->setItem(ii, 3, new QTableWidgetItem(tr("NULL")));
					ui.tableWidget->setItem(ii, 4, new QTableWidgetItem(tr("NULL")));
					QFont nullFont;
					QBrush nullColor(Qt::black);
					ui.tableWidget->item(ii, 2)->setFont(nullFont);
					ui.tableWidget->item(ii, 2)->setForeground(nullColor);
					ui.tableWidget->item(ii, 3)->setFont(nullFont);
					ui.tableWidget->item(ii, 3)->setForeground(nullColor);
					ui.tableWidget->item(ii, 4)->setFont(nullFont);
					ui.tableWidget->item(ii, 4)->setForeground(nullColor);
				}
				else
				{
					if (vec_temp.size() % 2 == 0)
					{
						median_err = (vec_temp[vec_temp.size() / 2] + vec_temp[vec_temp.size() / 2 - 1]) / 2.0;
					}
					else
					{
						median_err = vec_temp[(vec_temp.size() - 1) / 2];
					}
					mean_err /= (double)err_number;
					ui.tableWidget->setItem(ii, 2, new QTableWidgetItem(QString::number(mean_err, 'f', 4)));
					ui.tableWidget->setItem(ii, 3, new QTableWidgetItem(QString::number(median_err, 'f', 4)));
					ui.tableWidget->setItem(ii, 4, new QTableWidgetItem(QString::number(max_err, 'f', 4)));
					if (ui.reproject_error_limit_doubleSpinBox->isEnabled()
						&& mean_err >= ui.reproject_error_limit_doubleSpinBox->value())
					{
						QFont nullFont;
						nullFont.setBold(true);
						QBrush nullColor(Qt::red);
						ui.tableWidget->item(ii, 2)->setFont(nullFont);
						ui.tableWidget->item(ii, 2)->setForeground(nullColor);

					}
					else
					{
						QFont nullFont;
						QBrush nullColor(Qt::black);
						ui.tableWidget->item(ii, 2)->setFont(nullFont);
						ui.tableWidget->item(ii, 2)->setForeground(nullColor);
					}
					if (ui.reproject_error_limit_doubleSpinBox->isEnabled()
						&& max_err >= ui.reproject_error_limit_doubleSpinBox->value())
					{
						QFont nullFont;
						nullFont.setBold(true);
						QBrush nullColor(Qt::red);
						ui.tableWidget->item(ii, 4)->setFont(nullFont);
						ui.tableWidget->item(ii, 4)->setForeground(nullColor);
					}
					else
					{
						QFont nullFont;
						QBrush nullColor(Qt::black);
						ui.tableWidget->item(ii, 4)->setFont(nullFont);
						ui.tableWidget->item(ii, 4)->setForeground(nullColor);
					}
					if (ui.reproject_error_limit_doubleSpinBox->isEnabled()
						&& median_err >= ui.reproject_error_limit_doubleSpinBox->value())
					{
						QFont nullFont;
						nullFont.setBold(true);
						QBrush nullColor(Qt::red);
						ui.tableWidget->item(ii, 3)->setFont(nullFont);
						ui.tableWidget->item(ii, 3)->setForeground(nullColor);
					}
					else
					{
						QFont nullFont;
						QBrush nullColor(Qt::black);
						ui.tableWidget->item(ii, 3)->setFont(nullFont);
						ui.tableWidget->item(ii, 3)->setForeground(nullColor);
					}
				}
			}
			ui.tableWidget->setItem(ii, 5, new QTableWidgetItem(Image_serial_name[ii]));
			int cur_index = 255.0 / 2.0;
			cur_index = (cur_index < 0) ? 0 : cur_index;
			cur_index = (cur_index > 255) ? 255 : cur_index;
			auto back_color = QBrush(QColor(LOOKUPTABLE_JET[cur_index * 3] * 255,
				LOOKUPTABLE_JET[cur_index * 3 + 1] * 255, LOOKUPTABLE_JET[cur_index * 3 + 2] * 255, 80));
			ui.tableWidget->item(ii, 0)->setBackground(back_color);
			ui.tableWidget->item(ii, 1)->setBackground(back_color);
			ui.tableWidget->item(ii, 2)->setBackground(back_color);
			ui.tableWidget->item(ii, 3)->setBackground(back_color);
			ui.tableWidget->item(ii, 4)->setBackground(back_color);
			ui.tableWidget->item(ii, 5)->setBackground(back_color);
		}
	}    
	if (pri_vec_all.size() != 0)
	{
		double median_err;
		double max_err;
		double mean_err=0;
		int number_err = 0;
		while (pri_vec_all.size())
		{
			vec_all_temp.push_back(pri_vec_all.top());
			mean_err += pri_vec_all.top();
			number_err++;
			pri_vec_all.pop();
		}
		mean_err /= (double)number_err;
		if (vec_all_temp.size() % 2 == 0)
		{
			median_err = (vec_all_temp[vec_all_temp.size() / 2] + vec_all_temp[vec_all_temp.size() / 2 - 1]) / 2.0;
		}
		else
		{
			median_err = vec_all_temp[(vec_all_temp.size() - 1) / 2];
		}
		max_err = vec_all_temp[0];
		ui.summary_label->setText(tr("Mean:") + QString::number(mean_err, 'f', 4) + tr("; Median:") + QString::number(median_err, 'f', 4)
			+ tr("; Max:") + QString::number(max_err, 'f', 4));
	}
	else
	{
		ui.summary_label->setText(tr(""));
	}
	QStringList image_header_h;
	QStringList image_header_v;
	ui.image_error_tableView->setRowCount(0);
	ui.image_error_tableView->setColumnCount(0);
	ui.image_error_tableView->setHorizontalHeaderLabels(image_header_h);
	ui.image_error_tableView->setVerticalHeaderLabels(image_header_v);
	if (Result_err_every.size() != 0 && KeyPoint_code_serial.size() == Result_err_every.size())
	{
		ui.image_error_tableView->setRowCount(KeyPoint_code_serial.size() + 2);
		ui.image_error_tableView->setColumnCount(Result_xyz_serial.size() + 2);
		for (int ii = 0; ii < KeyPoint_code_serial.size() + 2; ii++)
		{
			for (int jj = 0; jj < Result_xyz_serial.size() + 2; jj++)
			{
				if (ii < KeyPoint_code_serial.size() || jj < Result_xyz_serial.size())
				{
					ui.image_error_tableView->setItem(ii, jj, new QTableWidgetItem("--"));
					QFont nullFont;
					QBrush nullColor(Qt::black);
					ui.image_error_tableView->item(ii, jj)->setFont(nullFont);
					ui.image_error_tableView->item(ii, jj)->setForeground(nullColor);
				}
				else
				{
					ui.image_error_tableView->setItem(ii, jj, new QTableWidgetItem(""));
					QFont nullFont;
					QBrush nullColor(Qt::black);
					ui.image_error_tableView->item(ii, jj)->setFont(nullFont);
					ui.image_error_tableView->item(ii, jj)->setForeground(nullColor);
				}
			}
		}
		for (int ii = 0; ii < Result_err_every.size(); ii++)
		{
			image_header_v << tr("Image-") + QString::number(ii + 1);
		}
		image_header_v << tr("Number");
		image_header_v << tr("Mean");
		for (int ii = 0; ii < Result_xyz_serial.size(); ii++)
		{
			image_header_h << tr("Coded Circle-") + QString::number(Result_xyz_serial[ii].code);
		}
		image_header_h << tr("Number");
		image_header_h << tr("Mean");
		ui.image_error_tableView->setHorizontalHeaderLabels(image_header_h);
		ui.image_error_tableView->setVerticalHeaderLabels(image_header_v);
		for (int ii = 0; ii < Result_err_every.size(); ii++)
		{
			for (int jj = 0; jj < Result_err_every[ii].size(); jj++)
			{
				if (Result_err_every[ii][jj].x == std::numeric_limits<double>::max() || Result_err_every[ii][jj].y == std::numeric_limits<double>::max())
				{
					continue;
				}
				bool exist = false;
				int cur_index = 0;
				for (int pp = 0; pp < Result_xyz_serial.size(); pp++)
				{
					if (KeyPoint_code_serial[ii][jj].code_num == Result_xyz_serial[pp].code)
					{
						exist = true;
						cur_index = pp;
						break;
					}
				}
				if (exist)
				{
					if (ui.image_error_tableView->item(ii, cur_index) == NULL ||
						(ui.image_error_tableView->item(ii, cur_index)))
					{
						ui.image_error_tableView->setItem(ii, cur_index, new QTableWidgetItem(QString::number(Result_err_every[ii][jj].x, 'f', 3)
							+ tr("/") + QString::number(Result_err_every[ii][jj].y, 'f', 3)
							+ tr("/") + QString::number(sqrt(Result_err_every[ii][jj].x * Result_err_every[ii][jj].x + Result_err_every[ii][jj].y * Result_err_every[ii][jj].y), 'f', 3)));
						if (ui.reproject_error_limit_doubleSpinBox->isEnabled() &&
							(Result_err_every[ii][jj].x >= ui.reproject_error_limit_doubleSpinBox->value() || Result_err_every[ii][jj].y >= ui.reproject_error_limit_doubleSpinBox->value())
							|| sqrt(Result_err_every[ii][jj].x * Result_err_every[ii][jj].x + Result_err_every[ii][jj].y * Result_err_every[ii][jj].y) >= ui.reproject_error_limit_doubleSpinBox->value())
						{
							QFont nullFont;
							nullFont.setBold(true);
							QBrush nullColor(Qt::red);
							ui.image_error_tableView->item(ii, cur_index)->setFont(nullFont);
							ui.image_error_tableView->item(ii, cur_index)->setForeground(nullColor);
						}
						else
						{
							QFont nullFont;
							QBrush nullColor(Qt::black);
							ui.image_error_tableView->item(ii, cur_index)->setFont(nullFont);
							ui.image_error_tableView->item(ii, cur_index)->setForeground(nullColor);
						}
					}
				}
			}
		}
		for (int pp = 0; pp < Result_xyz_serial.size(); pp++)
		{
			int index_number=0;
			double mean_x = 0, mean_y = 0, mean_r = 0;
			for (int ii = 0; ii < Result_err_every.size(); ii++)
			{
				for (int jj = 0; jj < Result_err_every[ii].size(); jj++)
				{
					if (KeyPoint_code_serial[ii][jj].code_num == Result_xyz_serial[pp].code)
					{
						if (Result_err_every[ii][jj].x == std::numeric_limits<double>::max() || Result_err_every[ii][jj].y == std::numeric_limits<double>::max())
						{
							continue;
						}
						index_number++;
						mean_x += abs(Result_err_every[ii][jj].x);
						mean_y += abs(Result_err_every[ii][jj].y);
						mean_r += sqrt(Result_err_every[ii][jj].x * Result_err_every[ii][jj].x+ Result_err_every[ii][jj].y* Result_err_every[ii][jj].y);
					}
				}
			}
			ui.image_error_tableView->setItem(Result_err_every.size(), pp, new QTableWidgetItem(QString::number(index_number)));
			QFont nullFont_number;
			QBrush nullColor_number(Qt::black);
			ui.image_error_tableView->item(Result_err_every.size(), pp)->setFont(nullFont_number);
			ui.image_error_tableView->item(Result_err_every.size(), pp)->setForeground(nullColor_number);
			if (index_number != 0)
			{
				mean_x /= (double)index_number;
				mean_y /= (double)index_number;
				mean_r /= (double)index_number;
				ui.image_error_tableView->setItem(Result_err_every.size() + 1, pp, new QTableWidgetItem(QString::number(mean_x, 'f', 3)
					+ tr("/") + QString::number(mean_y, 'f', 3)
					+ tr("/") + QString::number(mean_r, 'f', 3)));
				if (ui.reproject_error_limit_doubleSpinBox->isEnabled() &&
					(mean_x >= ui.reproject_error_limit_doubleSpinBox->value() || mean_y >= ui.reproject_error_limit_doubleSpinBox->value())
					|| mean_r >= ui.reproject_error_limit_doubleSpinBox->value())
				{
					QFont nullFont;
					nullFont.setBold(true);
					QBrush nullColor(Qt::red);
					ui.image_error_tableView->item(Result_err_every.size() + 1, pp)->setFont(nullFont);
					ui.image_error_tableView->item(Result_err_every.size() + 1, pp)->setForeground(nullColor);
				}
				else
				{
					QFont nullFont;
					QBrush nullColor(Qt::black);
					ui.image_error_tableView->item(Result_err_every.size() + 1, pp)->setFont(nullFont);
					ui.image_error_tableView->item(Result_err_every.size() + 1, pp)->setForeground(nullColor);
				}
			}
			else
			{
				ui.image_error_tableView->setItem(Result_err_every.size() + 1, pp, new QTableWidgetItem(tr("NULL")));
				QFont nullFont;
				QBrush nullColor(Qt::black);
				ui.image_error_tableView->item(Result_err_every.size() + 1, pp)->setFont(nullFont);
				ui.image_error_tableView->item(Result_err_every.size() + 1, pp)->setForeground(nullColor);
			}
		}
		for (int ii = 0; ii < Result_err_every.size(); ii++)
		{
			int index_number = 0;
			double mean_x = 0, mean_y = 0, mean_r = 0;
			for (int jj = 0; jj < Result_err_every[ii].size(); jj++)
			{
				if (Result_err_every[ii][jj].x == std::numeric_limits<double>::max() || Result_err_every[ii][jj].y == std::numeric_limits<double>::max())
				{
					continue;
				}
				index_number++;
				mean_x += abs(Result_err_every[ii][jj].x);
				mean_y += abs(Result_err_every[ii][jj].y);
				mean_r += sqrt(Result_err_every[ii][jj].x * Result_err_every[ii][jj].x + Result_err_every[ii][jj].y * Result_err_every[ii][jj].y);
			}
			ui.image_error_tableView->setItem(ii, Result_xyz_serial.size(), new QTableWidgetItem(QString::number(index_number)));
			QFont nullFont_number;
			QBrush nullColor_number(Qt::black);
			ui.image_error_tableView->item(ii, Result_xyz_serial.size())->setFont(nullFont_number);
			ui.image_error_tableView->item(ii, Result_xyz_serial.size())->setForeground(nullColor_number);
			if (index_number != 0)
			{
				mean_x /= (double)index_number;
				mean_y /= (double)index_number;
				mean_r /= (double)index_number;
				ui.image_error_tableView->setItem(ii, Result_xyz_serial.size() + 1, new QTableWidgetItem(QString::number(mean_x, 'f', 3)
					+ tr("/") + QString::number(mean_y, 'f', 3)
					+ tr("/") + QString::number(mean_r, 'f', 3)));
				if (ui.reproject_error_limit_doubleSpinBox->isEnabled() &&
					(mean_x >= ui.reproject_error_limit_doubleSpinBox->value() || mean_y >= ui.reproject_error_limit_doubleSpinBox->value())
					|| mean_r >= ui.reproject_error_limit_doubleSpinBox->value())
				{
					QFont nullFont;
					nullFont.setBold(true);
					QBrush nullColor(Qt::red);
					ui.image_error_tableView->item(ii, Result_xyz_serial.size() + 1)->setFont(nullFont);
					ui.image_error_tableView->item(ii, Result_xyz_serial.size() + 1)->setForeground(nullColor);
				}
				else
				{
					QFont nullFont;
					QBrush nullColor(Qt::black);
					ui.image_error_tableView->item(ii, Result_xyz_serial.size() + 1)->setFont(nullFont);
					ui.image_error_tableView->item(ii, Result_xyz_serial.size() + 1)->setForeground(nullColor);
				}
			}
			else
			{
				ui.image_error_tableView->setItem(ii, Result_xyz_serial.size() + 1, new QTableWidgetItem(tr("NULL")));
				QFont nullFont;
				QBrush nullColor(Qt::black);
				ui.image_error_tableView->item(ii, Result_xyz_serial.size() + 1)->setFont(nullFont);
				ui.image_error_tableView->item(ii, Result_xyz_serial.size() + 1)->setForeground(nullColor);
			}
		}
		
	}
	ui.image_error_tableView->resizeColumnsToContents();
	in_setting_keypoint_table = false;
}
void SFM_module::update_show_view_once(cv::Mat image, int show_index)
{
	cv::Mat cur_img = image;

	if (show_index < KeyPoint_code_serial.size())
	{
		if (!cur_img.empty())
		{
			double draw_R = std::numeric_limits<double>::max();
			double draw_R_min = std::numeric_limits<double>::max();
			if (KeyPoint_code_serial[show_index].size() == 1)
			{
				draw_R = sqrt(pow(KeyPoint_code_serial[show_index][0].r_a, 2)
					+ pow(KeyPoint_code_serial[show_index][0].r_b, 2));
			}
			else
			{
				for (int pp = 0; pp < KeyPoint_code_serial[show_index].size(); pp++)
				{
					for (int qq = pp + 1; qq < KeyPoint_code_serial[show_index].size(); qq++)
					{
						double dis_now = sqrt(pow(KeyPoint_code_serial[show_index][pp].x - KeyPoint_code_serial[show_index][qq].x, 2)
							+ pow(KeyPoint_code_serial[show_index][pp].y - KeyPoint_code_serial[show_index][qq].y, 2));
						draw_R = draw_R > dis_now ? dis_now : draw_R;
						double r_now = sqrt(pow(KeyPoint_code_serial[show_index][pp].r_a, 2)
							+ pow(KeyPoint_code_serial[show_index][pp].r_b, 2));
						draw_R_min = draw_R_min > r_now ? r_now : draw_R_min;
					}
				}
				draw_R /= 5;
				draw_R = draw_R < 10 ? 10 : draw_R;
				draw_R = draw_R > draw_R_min ? draw_R_min : draw_R;
			}
			double draw_thick = draw_R / 2;
			draw_thick = draw_thick < 1 ? 1 : draw_thick;
			if (Result_contours_serial.size() > show_index)
			{
				cv::drawContours(cur_img, Result_contours_serial[show_index], -1, cv::Scalar(0, 255, 0), draw_thick);
			}
			for (int pp = 0; pp < KeyPoint_code_serial[show_index].size(); pp++)
			{
				cv::putText(cur_img, std::to_string(KeyPoint_code_serial[show_index][pp].code_num)
					, cv::Point2f(KeyPoint_code_serial[show_index][pp].x, KeyPoint_code_serial[show_index][pp].y)
					, cv::FONT_HERSHEY_SIMPLEX, sqrt((double)cur_img.cols* (double)cur_img.rows/1000.0/1000.0), CV_RGB(255, 0, 0), 8);
			}
			ui.Detect_view->cam_ratio = (double)cur_img.cols / (double)cur_img.rows;
			ui.Detect_view->set_image(QImage((const unsigned char*)cur_img.data, cur_img.cols,
				cur_img.rows, cur_img.step, QImage::Format_BGR888).copy());
			ui.Detect_view->update();
		}
	}
}

void SFM_module::update_points_view()
{
	if (Result_xyz_serial.size() != 0)
	{
		ui.draw_setting_pushButton->setEnabled(true);
		in_setting_levelling_table = true;
		ui.levelling_tableWidget->setRowCount(0);
		ui.levelling_tableWidget->setRowCount(xyz_need_levelling_code_serial.size());
		for (int ii = 0; ii < xyz_need_levelling_code_serial.size(); ii++)
		{
			if (ui.levelling_tableWidget->item(ii, 0) == NULL ||
				(ui.levelling_tableWidget->item(ii, 0) &&
					ui.levelling_tableWidget->item(ii, 0)->text() == tr("")))
			{
				ui.levelling_tableWidget->setItem(ii, 0,
					new QTableWidgetItem(QString::number(xyz_need_levelling_code_serial[ii])));
				if (xyz_need_levelling_serial[ii])
				{
					ui.levelling_tableWidget->item(ii, 0)->setCheckState(Qt::Checked);
				}
				else
				{
					ui.levelling_tableWidget->item(ii, 0)->setCheckState(Qt::Unchecked);
				}
			}
		}
		ui.image_error_tableView->resizeColumnsToContents();

		in_setting_levelling_table = false;

		std::vector<sfm_3d_group> Result_xyz_serial_scale;
		for (int ii = 0; ii < Result_xyz_serial.size(); ii++)
		{
			Result_xyz_serial_scale.push_back(sfm_3d_group(Result_xyz_serial[ii].code,
				Result_xyz_serial[ii].weight,
				Result_xyz_serial[ii].x * scale_for_sfm, Result_xyz_serial[ii].y * scale_for_sfm, Result_xyz_serial[ii].z * scale_for_sfm, Result_xyz_serial[ii].error, Result_xyz_serial[ii].R * scale_for_sfm));
		}

		std::vector<cv::Point3d> points_need_fit;
		for (int ii = 0; ii < xyz_need_levelling_serial.size(); ii++)
		{
			if (xyz_need_levelling_serial[ii])
			{
				for (int jj = 0; jj < Result_xyz_serial_scale.size(); jj++)
				{
					if (Result_xyz_serial_scale[jj].code == xyz_need_levelling_code_serial[ii])
					{
						points_need_fit.push_back(cv::Point3d(
							Result_xyz_serial_scale[jj].x, Result_xyz_serial_scale[jj].y, Result_xyz_serial_scale[jj].z));
					}
				}
			}
		}

		if (points_need_fit.size() < 3)
		{
			Levelling_R = Eigen::Matrix3d::Zero(3, 3);
			Levelling_R(0, 0) = 1;
			Levelling_R(1, 1) = 1;
			Levelling_R(2, 2) = 1;
			Levelling_T = Eigen::Vector3d(0, 0, 0);
		}
		else
		{
			int num_s = points_need_fit.size();
			Eigen::MatrixXd mleft = Eigen::MatrixXd::Zero(num_s, 3);
			double mean_x = 0, mean_y = 0, mean_z = 0;
			double mean_xx = 0, mean_yy = 0, mean_zz = 0;
			double mean_xy = 0, mean_yz = 0, mean_zx = 0;
			for (int i = 0; i < num_s; ++i)
			{
				mean_x += points_need_fit[i].x;
				mean_y += points_need_fit[i].y;
				mean_z += points_need_fit[i].z;
				mean_xx += points_need_fit[i].x * points_need_fit[i].x;
				mean_yy += points_need_fit[i].y * points_need_fit[i].y;
				mean_zz += points_need_fit[i].z * points_need_fit[i].z;
				mean_xy += points_need_fit[i].x * points_need_fit[i].y;
				mean_yz += points_need_fit[i].y * points_need_fit[i].z;
				mean_zx += points_need_fit[i].z * points_need_fit[i].x;
			}
			mean_x = mean_x / (double)num_s;
			mean_y = mean_y / (double)num_s;
			mean_xx = mean_xx / (double)num_s;
			mean_yy = mean_yy / (double)num_s;
			mean_xy = mean_xy / (double)num_s;
			mean_z = mean_z / (double)num_s;
			mean_zz = mean_zz / (double)num_s;
			mean_yz = mean_yz / (double)num_s;
			mean_zx = mean_zx / (double)num_s;
			Eigen::Matrix3d eA;
			eA(0, 0) = mean_xx - mean_x * mean_x; eA(0, 1) = mean_xy - mean_x * mean_y; eA(0, 2) = mean_zx - mean_x * mean_z;
			eA(1, 0) = mean_xy - mean_x * mean_y; eA(1, 1) = mean_yy - mean_y * mean_y; eA(1, 2) = mean_yz - mean_y * mean_z;
			eA(2, 0) = mean_zx - mean_x * mean_z; eA(2, 1) = mean_yz - mean_y * mean_z; eA(2, 2) = mean_zz - mean_z * mean_z;
			Eigen::EigenSolver<Eigen::Matrix3d> eMat(eA);
			Eigen::Matrix3d eD = eMat.pseudoEigenvalueMatrix();
			Eigen::Matrix3d eV = eMat.pseudoEigenvectors();

			/* the eigenvector corresponding to the minimum eigenvalue */
			double eD1 = eD(0, 0);
			double eD2 = eD(1, 1);
			double eD3 = eD(2, 2);
			int minNumber = 0;
			if ((abs(eD2) <= abs(eD1)) && (abs(eD2) <= abs(eD3)))
			{
				minNumber = 1;
			}
			if ((abs(eD3) <= abs(eD1)) && (abs(eD3) <= abs(eD2)))
			{
				minNumber = 2;
			}
			double A = eV(0, minNumber);
			double B = eV(1, minNumber);
			double C = eV(2, minNumber);
			double D = -(A * mean_x + B * mean_y + C * mean_z);

			/* result */
			if (ui.z_reversal_checkBox->isChecked())
			{
				A *= -1;
				B *= -1;
				C *= -1;
				D *= -1;
			}
			Eigen::Vector3d norm_p;
			norm_p.x() = A / sqrt(A * A + B * B + C * C);
			norm_p.y() = B / sqrt(A * A + B * B + C * C);
			norm_p.z() = C / sqrt(A * A + B * B + C * C);

			Eigen::Vector3d norm_xy = Eigen::Vector3d(0, 0, 1);

			Levelling_R = Eigen::Quaterniond::FromTwoVectors(norm_p, norm_xy).toRotationMatrix();
			double sum_x = 0, sum_y = 0, sum_z = 0;
			for (int i = 0; i < num_s; ++i)
			{
				Eigen::Vector3d point_old = Eigen::Vector3d(points_need_fit[i].x, points_need_fit[i].y, points_need_fit[i].z);
				Eigen::Vector3d point_now = Levelling_R * point_old;
				sum_x += point_now.x();
				sum_y += point_now.y();
				sum_z += point_now.z();
			}
			sum_x /= (double)num_s;
			sum_y /= (double)num_s;
			sum_z /= (double)num_s;
			Levelling_T = Eigen::Vector3d(-sum_x, -sum_y, -sum_z);
		}
		Result_xyz_after_op_serial.clear();
		for (std::size_t ii = 0; ii < Result_xyz_serial_scale.size(); ++ii)
		{
			Eigen::Vector3d point_old = Eigen::Vector3d(Result_xyz_serial_scale[ii].x, Result_xyz_serial_scale[ii].y, Result_xyz_serial_scale[ii].z);
			Eigen::Vector3d point_now = Levelling_R * point_old+ Levelling_T;
			Result_xyz_after_op_serial.push_back(sfm_3d_group(Result_xyz_serial_scale[ii].code,
				Result_xyz_serial_scale[ii].weight,
				point_now.x(), point_now.y(), point_now.z(), Result_xyz_serial_scale[ii].error, Result_xyz_serial_scale[ii].R));
		}

		double point_size = std::numeric_limits<double>::max();
		for (std::size_t ii = 0; ii < Result_xyz_after_op_serial.size(); ++ii)
		{
			for (std::size_t jj = ii + 1; jj < Result_xyz_after_op_serial.size(); ++jj)
			{
				double r = sqrt(pow((Result_xyz_after_op_serial[ii].x - Result_xyz_after_op_serial[jj].x), 2)
					+ pow((Result_xyz_after_op_serial[ii].y - Result_xyz_after_op_serial[jj].y), 2)
					+ pow((Result_xyz_after_op_serial[ii].z - Result_xyz_after_op_serial[jj].z), 2));
				point_size = point_size > r ? r : point_size;
			}
		}

		window_for_vtk->RemoveRenderer(viewer_for_vtk->getRendererCollection()->GetFirstRenderer());
		viewer_for_vtk.reset(new pcl::visualization::PCLVisualizer("viewer", false));
		window_for_vtk->AddRenderer(viewer_for_vtk->getRendererCollection()->GetFirstRenderer());
		ui.VTK_xyz_widget->setRenderWindow(window_for_vtk.Get());
		viewer_for_vtk->setupInteractor(ui.VTK_xyz_widget->interactor(), ui.VTK_xyz_widget->renderWindow());
		viewer_for_vtk->setBackgroundColor(1, 1, 1);
		ui.VTK_xyz_widget->update();

		cloud_for_vtk.reset(new pcl::PointCloud<pcl::PointXYZRGBA>);
		cloud_for_vtk->points.resize(Result_xyz_after_op_serial.size());
		for (std::size_t ii = 0; ii < cloud_for_vtk->points.size(); ++ii)
		{
			cloud_for_vtk->points[ii].x = Result_xyz_after_op_serial[ii].x;
			cloud_for_vtk->points[ii].y = Result_xyz_after_op_serial[ii].y;
			cloud_for_vtk->points[ii].z = Result_xyz_after_op_serial[ii].z;
		}
		Points_xyz_range range_now = VTK_plot::get_xyz_range(cloud_for_vtk);

		PCL_setting_window->in_setting = true;
		PCL_setting_window->ui.min_x_doubleSpinBox->setMinimum(range_now.x_min - (range_now.x_max - range_now.x_min) * 0.25);
		PCL_setting_window->ui.max_x_doubleSpinBox->setMinimum(range_now.x_min - (range_now.x_max - range_now.x_min) * 0.25);
		PCL_setting_window->ui.min_x_doubleSpinBox->setMaximum(range_now.x_max + (range_now.x_max - range_now.x_min) * 0.25);
		PCL_setting_window->ui.max_x_doubleSpinBox->setMaximum(range_now.x_max + (range_now.x_max - range_now.x_min) * 0.25);

		PCL_setting_window->ui.min_y_doubleSpinBox->setMinimum(range_now.y_min - (range_now.y_max - range_now.y_min) * 0.25);
		PCL_setting_window->ui.max_y_doubleSpinBox->setMinimum(range_now.y_min - (range_now.y_max - range_now.y_min) * 0.25);
		PCL_setting_window->ui.min_y_doubleSpinBox->setMaximum(range_now.y_max + (range_now.y_max - range_now.y_min) * 0.25);
		PCL_setting_window->ui.max_y_doubleSpinBox->setMaximum(range_now.y_max + (range_now.y_max - range_now.y_min) * 0.25);

		PCL_setting_window->ui.min_z_doubleSpinBox->setMinimum(range_now.z_min - (range_now.z_max - range_now.z_min) * 0.25);
		PCL_setting_window->ui.max_z_doubleSpinBox->setMinimum(range_now.z_min - (range_now.z_max - range_now.z_min) * 0.25);
		PCL_setting_window->ui.min_z_doubleSpinBox->setMaximum(range_now.z_max + (range_now.z_max - range_now.z_min) * 0.25);
		PCL_setting_window->ui.max_z_doubleSpinBox->setMaximum(range_now.z_max + (range_now.z_max - range_now.z_min) * 0.25);

		if (PCL_setting_window->ui.auto_x_checkbox->isChecked())
		{
			PCL_setting_window->ui.min_x_doubleSpinBox->setValue(range_now.x_min);
			PCL_setting_window->ui.max_x_doubleSpinBox->setValue(range_now.x_max);
		}
		if (PCL_setting_window->ui.auto_y_checkbox->isChecked())
		{
			PCL_setting_window->ui.min_y_doubleSpinBox->setValue(range_now.y_min);
			PCL_setting_window->ui.max_y_doubleSpinBox->setValue(range_now.y_max);
		}
		if (PCL_setting_window->ui.auto_z_checkbox->isChecked())
		{
			PCL_setting_window->ui.min_z_doubleSpinBox->setValue(range_now.z_min);
			PCL_setting_window->ui.max_z_doubleSpinBox->setValue(range_now.z_max);
		}

		PCL_setting_window->in_setting = false;
		pcl::PointCloud<pcl::PointXYZRGBA>::Ptr cloud_for_vtk_copy;
		cloud_for_vtk_copy.reset(new pcl::PointCloud<pcl::PointXYZRGBA>);
		std::vector<bool> show_points;
		for (std::size_t ii = 0; ii < cloud_for_vtk->points.size(); ++ii)
		{
			if ((PCL_setting_window->ui.auto_x_checkbox->isChecked() ||
				(cloud_for_vtk->points[ii].x > PCL_setting_window->ui.min_x_doubleSpinBox->value() && cloud_for_vtk->points[ii].x < PCL_setting_window->ui.max_x_doubleSpinBox->value()))&&
				(PCL_setting_window->ui.auto_y_checkbox->isChecked() ||
					(cloud_for_vtk->points[ii].y > PCL_setting_window->ui.min_y_doubleSpinBox->value() && cloud_for_vtk->points[ii].y < PCL_setting_window->ui.max_y_doubleSpinBox->value())) &&
				(PCL_setting_window->ui.auto_z_checkbox->isChecked() ||
					(cloud_for_vtk->points[ii].z > PCL_setting_window->ui.min_z_doubleSpinBox->value() && cloud_for_vtk->points[ii].z < PCL_setting_window->ui.max_z_doubleSpinBox->value())))
			{
				show_points.push_back(true);
				cloud_for_vtk_copy->points.push_back(pcl::PointXYZRGBA(cloud_for_vtk->points[ii].x, cloud_for_vtk->points[ii].y, cloud_for_vtk->points[ii].z));
			}
			else
			{
				show_points.push_back(false);
			}
		}

		range_now = VTK_plot::get_xyz_range(cloud_for_vtk_copy);

		PCL_setting_window->in_setting = true;
		switch (PCL_setting_window->ui.color_xyz_comboBox->currentIndex())
		{
		case 0:
			PCL_setting_window->ui.min_xyz_lut_doubleSpinBox->setMinimum(range_now.x_min - (range_now.x_max - range_now.x_min) * 0.25);
			PCL_setting_window->ui.max_xyz_lut_doubleSpinBox->setMinimum(range_now.x_min - (range_now.x_max - range_now.x_min) * 0.25);
			PCL_setting_window->ui.min_xyz_lut_doubleSpinBox->setMaximum(range_now.x_max + (range_now.x_max - range_now.x_min) * 0.25);
			PCL_setting_window->ui.max_xyz_lut_doubleSpinBox->setMaximum(range_now.x_max + (range_now.x_max - range_now.x_min) * 0.25);
			if (PCL_setting_window->ui.auto_lut_xyz_checkbox->isChecked())
			{
				PCL_setting_window->ui.min_xyz_lut_doubleSpinBox->setValue(range_now.x_min);
				PCL_setting_window->ui.max_xyz_lut_doubleSpinBox->setValue(range_now.x_max);
			}
			break;
		case 1:
			PCL_setting_window->ui.min_xyz_lut_doubleSpinBox->setMinimum(range_now.y_min - (range_now.y_max - range_now.y_min) * 0.25);
			PCL_setting_window->ui.max_xyz_lut_doubleSpinBox->setMinimum(range_now.y_min - (range_now.y_max - range_now.y_min) * 0.25);
			PCL_setting_window->ui.min_xyz_lut_doubleSpinBox->setMaximum(range_now.y_max + (range_now.y_max - range_now.y_min) * 0.25);
			PCL_setting_window->ui.max_xyz_lut_doubleSpinBox->setMaximum(range_now.y_max + (range_now.y_max - range_now.y_min) * 0.25);
			if (PCL_setting_window->ui.auto_lut_xyz_checkbox->isChecked())
			{
				PCL_setting_window->ui.min_xyz_lut_doubleSpinBox->setValue(range_now.z_min);
				PCL_setting_window->ui.max_xyz_lut_doubleSpinBox->setValue(range_now.z_max);
			}
			break;
		case 2:
			PCL_setting_window->ui.min_xyz_lut_doubleSpinBox->setMinimum(range_now.z_min - (range_now.z_max - range_now.z_min) * 0.25);
			PCL_setting_window->ui.max_xyz_lut_doubleSpinBox->setMinimum(range_now.z_min - (range_now.z_max - range_now.z_min) * 0.25);
			PCL_setting_window->ui.min_xyz_lut_doubleSpinBox->setMaximum(range_now.z_max + (range_now.z_max - range_now.z_min) * 0.25);
			PCL_setting_window->ui.max_xyz_lut_doubleSpinBox->setMaximum(range_now.z_max + (range_now.z_max - range_now.z_min) * 0.25);
			if (PCL_setting_window->ui.auto_lut_xyz_checkbox->isChecked())
			{
				PCL_setting_window->ui.min_xyz_lut_doubleSpinBox->setValue(range_now.z_min);
				PCL_setting_window->ui.max_xyz_lut_doubleSpinBox->setValue(range_now.z_max);
			}
			break;
		case 3:
		{
			std::vector<double> weight_vec;
			for (std::size_t ii = 0; ii < Result_xyz_after_op_serial.size(); ++ii)
			{
				weight_vec.push_back(Result_xyz_after_op_serial[ii].weight);
			}
			Points_xyz_range range_weight = VTK_plot::get_vec_range(weight_vec);
			PCL_setting_window->ui.min_xyz_lut_doubleSpinBox->setMinimum(range_weight.x_min);
			PCL_setting_window->ui.max_xyz_lut_doubleSpinBox->setMinimum(range_weight.x_min);
			PCL_setting_window->ui.min_xyz_lut_doubleSpinBox->setMaximum(range_weight.x_max);
			PCL_setting_window->ui.max_xyz_lut_doubleSpinBox->setMaximum(range_weight.x_max);
			if (PCL_setting_window->ui.auto_lut_xyz_checkbox->isChecked())
			{
				PCL_setting_window->ui.min_xyz_lut_doubleSpinBox->setValue(range_weight.x_min);
				PCL_setting_window->ui.max_xyz_lut_doubleSpinBox->setValue(range_weight.x_max);
			}
			break;
		}
		default:
			break;
		}

		PCL_setting_window->in_setting = false;

		point_size = point_size > ((range_now.x_max - range_now.x_min) / 5.0) ? ((range_now.x_max - range_now.x_min) / 5.0) : point_size;
		point_size = point_size > ((range_now.y_max - range_now.y_min) / 5.0) ? ((range_now.y_max - range_now.y_min) / 5.0) : point_size;
		point_size = point_size > ((range_now.z_max - range_now.z_min) / 5.0) ? ((range_now.z_max - range_now.z_min) / 5.0) : point_size;
		
		VTK_plot::set_grid_mode(viewer_for_vtk, cubeaxes_for_vtk, range_now, PCL_setting_window->ui.grid_on_checkbox->isChecked());

		Points_xyz_range range_now_copy = range_now;
		if (!PCL_setting_window->ui.auto_lut_xyz_checkbox->isChecked())
		{
			range_now_copy.x_min = PCL_setting_window->ui.min_xyz_lut_doubleSpinBox->value();
			range_now_copy.y_min = PCL_setting_window->ui.min_xyz_lut_doubleSpinBox->value();
			range_now_copy.z_min = PCL_setting_window->ui.min_xyz_lut_doubleSpinBox->value();
			range_now_copy.x_max = PCL_setting_window->ui.max_xyz_lut_doubleSpinBox->value();
			range_now_copy.y_max = PCL_setting_window->ui.max_xyz_lut_doubleSpinBox->value();
			range_now_copy.z_max = PCL_setting_window->ui.max_xyz_lut_doubleSpinBox->value();
		}
		VTK_plot::set_lut_mode(lut_for_vtk, PCL_setting_window->ui.lookuptable_color_comboBox->currentIndex());
		if (PCL_setting_window->ui.color_xyz_comboBox->currentIndex() < 3)
		{
			VTK_plot::set_clouds_color_xyz(cloud_for_vtk_copy, lut_for_vtk, range_now_copy, PCL_setting_window->ui.color_xyz_comboBox->currentIndex());

			if (PCL_setting_window->ui.colorbar_on_checkbox->isChecked() && PCL_setting_window->ui.color_xyz_comboBox->currentIndex() == 0)
			{
				VTK_plot::set_lut_show(viewer_for_vtk, lut_for_vtk, scalarBarActor_for_vtk, range_now_copy.x_min, range_now_copy.x_max, tr("X").toLocal8Bit().toStdString().c_str());
			}
			if (PCL_setting_window->ui.colorbar_on_checkbox->isChecked() && PCL_setting_window->ui.color_xyz_comboBox->currentIndex() == 1)
			{
				VTK_plot::set_lut_show(viewer_for_vtk, lut_for_vtk, scalarBarActor_for_vtk, range_now_copy.y_min, range_now_copy.y_max, tr("Y").toLocal8Bit().toStdString().c_str());
			}
			if (PCL_setting_window->ui.colorbar_on_checkbox->isChecked() && PCL_setting_window->ui.color_xyz_comboBox->currentIndex() == 2)
			{
				VTK_plot::set_lut_show(viewer_for_vtk, lut_for_vtk, scalarBarActor_for_vtk, range_now_copy.z_min, range_now_copy.z_max, tr("Z").toLocal8Bit().toStdString().c_str());
			}
		}
		else if (PCL_setting_window->ui.color_xyz_comboBox->currentIndex() == 3)
		{
			std::vector<double> weight_vec;
			for (std::size_t ii = 0; ii < Result_xyz_after_op_serial.size(); ++ii)
			{
				weight_vec.push_back(Result_xyz_after_op_serial[ii].weight);
			}
			Points_xyz_range range_weight = VTK_plot::get_vec_range(weight_vec);
			VTK_plot::set_clouds_color_other(cloud_for_vtk_copy, weight_vec, lut_for_vtk, range_weight);

			if (PCL_setting_window->ui.colorbar_on_checkbox->isChecked())
			{
				VTK_plot::set_lut_show(viewer_for_vtk, lut_for_vtk, scalarBarActor_for_vtk, range_weight.x_min, range_weight.x_max, tr("View Count").toLocal8Bit().toStdString().c_str());
			}
		}

		viewer_for_vtk->getRendererCollection()->GetFirstRenderer()->GetActors();
		viewer_for_vtk->getRendererCollection()->GetFirstRenderer()->ResetCamera();
		viewer_for_vtk->getRendererCollection()->GetFirstRenderer()->ResetCameraClippingRange();
		viewer_for_vtk->removeAllPointClouds();
		viewer_for_vtk->removeAllShapes();
		viewer_for_vtk->addPointCloud(cloud_for_vtk_copy, "cloud_xyz");
		viewer_for_vtk->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, PCL_setting_window->ui.point_size_spin->value(), "cloud_xyz");
		viewer_for_vtk->spinOnce();
;

		viewer_for_vtk->updatePointCloud(cloud_for_vtk_copy, "cloud_xyz");

		if (Result_R.size() == Image_serial_name.size() && Result_T.size() == Image_serial_name.size())
		{
			Result_R_after_op.clear();
			Result_T_after_op.clear();
			for (int ii = 0; ii < Result_xyz_after_op_serial.size(); ii++)
			{
				if (ii < Result_T.size() && ii < Result_R.size())
				{
					Eigen::Vector3d point_old = Eigen::Vector3d(Result_T[ii].x() * scale_for_sfm, Result_T[ii].y() * scale_for_sfm, Result_T[ii].z() * scale_for_sfm);
					Eigen::Vector3d point_now = Levelling_R * point_old + Levelling_T;
					Result_T_after_op.push_back(point_now);
					Eigen::Matrix3d r_now = Levelling_R * Result_R[ii].matrix();
					Result_R_after_op.push_back(Eigen::Quaterniond(r_now));
					Result_R_after_op[Result_R_after_op.size() - 1].normalize();
				}

			}
		}
		if (Result_T.size() == Image_serial_name.size() && PCL_setting_window->ui.path_checkBox->isChecked())
		{
			pcl::PointXYZ point_begin, point_end;
			bool find_first = false;
			int number_now = 0;
			for (int ii = 0; ii < Image_serial_name.size(); ii++)
			{
				if (!(Result_T[ii].x() == 0 && Result_T[ii].y() == 0 && Result_T[ii].z() == 0))
				{
					if (!find_first)
					{
						Eigen::Vector3d point_old = Eigen::Vector3d(Result_T[ii].x() * scale_for_sfm, Result_T[ii].y() * scale_for_sfm, Result_T[ii].z() * scale_for_sfm);
						Eigen::Vector3d point_now = Levelling_R * point_old + Levelling_T;
						point_begin.x = point_now.x();
						point_begin.y = point_now.y();
						point_begin.z = point_now.z();
						find_first = true;
					}
					else
					{
						Eigen::Vector3d point_old = Eigen::Vector3d(Result_T[ii].x() * scale_for_sfm, Result_T[ii].y() * scale_for_sfm, Result_T[ii].z() * scale_for_sfm);
						Eigen::Vector3d point_now = Levelling_R * point_old + Levelling_T;
						point_end.x = point_now.x();
						point_end.y = point_now.y();
						point_end.z = point_now.z();
						number_now++;
						QString ar_name = "line_" + QString::number(number_now);
						viewer_for_vtk->addLine<pcl::PointXYZ>(point_begin, point_end, 0, 0, 0, ar_name.toLocal8Bit().toStdString().c_str());
						viewer_for_vtk->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, PCL_setting_window->ui.point_size_spin->value() / 2.0, ar_name.toStdString().c_str());
						point_begin = point_end;
					}
				}
			}
		}
		if (Result_Ori.size() != 0 && Result_xyz_after_op_serial.size() == Result_Ori.size())
		{
			Result_Ori_after_op.clear();
			for (int ii = 0; ii < Result_xyz_after_op_serial.size(); ii++)
			{
				Eigen::Vector3d point_old =
					Eigen::Vector3d(Result_xyz_serial_scale[ii].x + Result_Ori[ii].x() * scale_for_sfm,
						Result_xyz_serial_scale[ii].y + Result_Ori[ii].y() * scale_for_sfm,
						Result_xyz_serial_scale[ii].z + Result_Ori[ii].z() * scale_for_sfm);
				Eigen::Vector3d point_now = Levelling_R * point_old + Levelling_T;
				pcl::PointXYZ point_arrow;
				point_arrow.x = Result_xyz_after_op_serial[ii].x;
				point_arrow.y = Result_xyz_after_op_serial[ii].y;
				point_arrow.z = Result_xyz_after_op_serial[ii].z;
				pcl::PointXYZ point_arrow_end;

				Eigen::Vector3d vec_now = Eigen::Vector3d(point_now.x() - Result_xyz_after_op_serial[ii].x,
					point_now.y() - Result_xyz_after_op_serial[ii].y,
					point_now.z() - Result_xyz_after_op_serial[ii].z);
				vec_now = vec_now / sqrt(vec_now.x() * vec_now.x() + vec_now.y() * vec_now.y() + vec_now.z() * vec_now.z());
				Result_Ori_after_op.push_back(vec_now);
			}
		}
		if (Result_Ori.size() !=0 && Result_xyz_after_op_serial.size() == Result_Ori.size() && PCL_setting_window->ui.ori_checkBox->isChecked())
		{
			pcl::PointXYZ point_begin, point_end;
			int number_arrow_now = 0;
			for (int ii = 0; ii < Result_xyz_after_op_serial.size(); ii++)
			{
				if (show_points[ii])
				{
					Eigen::Vector3d point_old = 
						Eigen::Vector3d(Result_xyz_serial_scale[ii].x + Result_Ori[ii].x() * scale_for_sfm,
							Result_xyz_serial_scale[ii].y + Result_Ori[ii].y() * scale_for_sfm,
							Result_xyz_serial_scale[ii].z + Result_Ori[ii].z() * scale_for_sfm);
					Eigen::Vector3d point_now = Levelling_R * point_old + Levelling_T;
					pcl::PointXYZ point_arrow;
					point_arrow.x =Result_xyz_after_op_serial[ii].x;
					point_arrow.y = Result_xyz_after_op_serial[ii].y;
					point_arrow.z = Result_xyz_after_op_serial[ii].z;
					pcl::PointXYZ point_arrow_end;



					point_arrow_end.x = point_arrow.x + Result_Ori_after_op[ii].x()* PCL_setting_window->ui.ori_scale_spin->value();
					point_arrow_end.y = point_arrow.y + Result_Ori_after_op[ii].y()* PCL_setting_window->ui.ori_scale_spin->value();
					point_arrow_end.z = point_arrow.z + Result_Ori_after_op[ii].z()* PCL_setting_window->ui.ori_scale_spin->value();
					QString ar_name = "arrowline_" + QString::number(number_arrow_now);
					viewer_for_vtk->addLine<pcl::PointXYZ>(point_arrow, point_arrow_end, 0, 0, 0, ar_name.toLocal8Bit().toStdString().c_str());
					viewer_for_vtk->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, PCL_setting_window->ui.point_size_spin->value() / 2.0, ar_name.toLocal8Bit().toStdString().c_str());
					number_arrow_now++;
				}
			}
		}


		if (PCL_setting_window->ui.code_checkBox->isChecked())
		{
			pcl::PointXYZ point_pos;
			int number_now = 0;
			for (int ii = 0; ii < Result_xyz_after_op_serial.size(); ii++)
			{
				if (show_points[ii])
				{
					point_pos.x = Result_xyz_after_op_serial[ii].x;
					point_pos.y = Result_xyz_after_op_serial[ii].y;
					point_pos.z = Result_xyz_after_op_serial[ii].z;
					QString ar_name = "text_" + QString::number(number_now);
					viewer_for_vtk->addText3D(std::to_string(Result_xyz_after_op_serial[ii].code), point_pos, PCL_setting_window->ui.text_scale_spin->value(), 0, 0, 0, ar_name.toLocal8Bit().toStdString().c_str());
					number_now++;
				}
			}
		}

		viewer_for_vtk->initCameraParameters();
		viewer_for_vtk->resetCamera();
		window_for_vtk->Render();
		ui.VTK_xyz_widget->renderWindow()->Render();
		viewer_for_vtk->setBackgroundColor(1, 1, 1);
		ui.VTK_xyz_widget->update();
		ui.VTK_xyz_widget->show();

	}
	else
	{
		window_for_vtk->RemoveRenderer(viewer_for_vtk->getRendererCollection()->GetFirstRenderer());
		viewer_for_vtk.reset(new pcl::visualization::PCLVisualizer("viewer", false));
		window_for_vtk->AddRenderer(viewer_for_vtk->getRendererCollection()->GetFirstRenderer());
		ui.VTK_xyz_widget->setRenderWindow(window_for_vtk.Get());
		viewer_for_vtk->setupInteractor(ui.VTK_xyz_widget->interactor(), ui.VTK_xyz_widget->renderWindow());
		viewer_for_vtk->setBackgroundColor(1, 1, 1);
		viewer_for_vtk->initCameraParameters();
		viewer_for_vtk->resetCamera();
		window_for_vtk->Render();
		ui.VTK_xyz_widget->renderWindow()->Render();
		viewer_for_vtk->setBackgroundColor(1, 1, 1);
		ui.VTK_xyz_widget->update();
		ui.VTK_xyz_widget->show();
	}
}
void SFM_module::update_show_view()
{
	if (in_calculation)
	{
		return;
	}
	if (Image_serial_name.size() == 0)
	{
		ui.show_detect_spinBox->setEnabled(false);
		ui.show_detect_horizontalSlider->setEnabled(false);
		ui.show_detect_spinBox->setMinimum(1);
		ui.show_detect_spinBox->setMaximum(1);
		ui.show_detect_horizontalSlider->setMinimum(1);
		ui.show_detect_horizontalSlider->setMaximum(1);
		ui.Detect_view->clear_image();
		ui.Detect_view->update();
		return;
	}
	else
	{
		ui.show_detect_spinBox->setEnabled(true);
		ui.show_detect_horizontalSlider->setEnabled(true);
		ui.show_detect_spinBox->setMinimum(1);
		ui.show_detect_spinBox->setMaximum(Image_serial_name.size());
		ui.show_detect_horizontalSlider->setMinimum(1);
		ui.show_detect_horizontalSlider->setMaximum(Image_serial_name.size());
	}
	if (Image_serial_name.size() == 1)
	{
		ui.show_detect_spinBox->setEnabled(false);
		ui.show_detect_horizontalSlider->setEnabled(false);
	}
	else
	{
		ui.show_detect_spinBox->setEnabled(true);
		ui.show_detect_horizontalSlider->setEnabled(true);
	}
	auto a = Image_serial_name[ui.show_detect_spinBox->value() - 1].toLocal8Bit().data();
	cv::Mat cur_img = cv::imread(Image_serial_name[ui.show_detect_spinBox->value() - 1].toLocal8Bit().data());

	if ((ui.show_detect_spinBox->value() - 1) < KeyPoint_code_serial.size())
	{
		if (cur_img.empty())
		{
			if (ImageWidth_serial[ui.show_detect_spinBox->value() - 1] > 0 &&
				ImageHeight_serial[ui.show_detect_spinBox->value() - 1] > 0)
			{
				cur_img = cv::Mat(ImageHeight_serial[ui.show_detect_spinBox->value() - 1]
					, ImageWidth_serial[ui.show_detect_spinBox->value() - 1], CV_8UC3, cv::Scalar(255, 255, 255));
			}
			else
			{
				ui.Detect_view->clear_image();
				ui.Detect_view->update();
				return;
			}
		}
		double draw_R = std::numeric_limits<double>::max();
		double draw_R_min = std::numeric_limits<double>::max();
		if (KeyPoint_code_serial[ui.show_detect_spinBox->value() - 1].size() == 1)
		{
			draw_R = sqrt(pow(KeyPoint_code_serial[ui.show_detect_spinBox->value() - 1][0].r_a, 2)
				+ pow(KeyPoint_code_serial[ui.show_detect_spinBox->value() - 1][0].r_b, 2));
		}
		else
		{
			for (int pp = 0; pp < KeyPoint_code_serial[ui.show_detect_spinBox->value() - 1].size(); pp++)
			{
				for (int qq = pp + 1; qq < KeyPoint_code_serial[ui.show_detect_spinBox->value() - 1].size(); qq++)
				{
					double dis_now = sqrt(pow(KeyPoint_code_serial[ui.show_detect_spinBox->value() - 1][pp].x - KeyPoint_code_serial[ui.show_detect_spinBox->value() - 1][qq].x, 2)
						+ pow(KeyPoint_code_serial[ui.show_detect_spinBox->value() - 1][pp].y - KeyPoint_code_serial[ui.show_detect_spinBox->value() - 1][qq].y, 2));
					draw_R = draw_R > dis_now ? dis_now : draw_R;
					double r_now = sqrt(pow(KeyPoint_code_serial[ui.show_detect_spinBox->value() - 1][pp].r_a, 2)
						+ pow(KeyPoint_code_serial[ui.show_detect_spinBox->value() - 1][pp].r_b, 2));
					draw_R_min = draw_R_min > r_now ? r_now : draw_R_min;
				}
			}
			draw_R /= 5;
			draw_R = draw_R < 10 ? 10 : draw_R;
			draw_R = draw_R > draw_R_min ? draw_R_min : draw_R;
		}
		double draw_thick = draw_R / 2;
		draw_thick = draw_thick < 1 ? 1 : draw_thick;
		if (Result_contours_serial.size() > (ui.show_detect_spinBox->value() - 1))
		{
			cv::drawContours(cur_img, Result_contours_serial[ui.show_detect_spinBox->value() - 1], -1, cv::Scalar(0, 255, 0), draw_thick);
		}
		for (int pp = 0; pp < KeyPoint_code_serial[ui.show_detect_spinBox->value() - 1].size(); pp++)
		{/*
			circle(cur_img, cv::Point2f(KeyPoint_code_serial[ui.show_detect_spinBox->value() - 1][pp].x, KeyPoint_code_serial[ui.show_detect_spinBox->value() - 1][pp].y)
				, draw_R, cv::Scalar(0, 255, 0), draw_thick);*/
			cv::putText(cur_img, std::to_string(KeyPoint_code_serial[ui.show_detect_spinBox->value() - 1][pp].code_num)
				, cv::Point2f(KeyPoint_code_serial[ui.show_detect_spinBox->value() - 1][pp].x, KeyPoint_code_serial[ui.show_detect_spinBox->value() - 1][pp].y)
				, cv::FONT_HERSHEY_SIMPLEX, sqrt((double)cur_img.cols * (double)cur_img.rows / 1000.0 / 1000.0), CV_RGB(255, 0, 0), 8);
		}
	}
	else
	{
		if (cur_img.empty())
		{
			ui.Detect_view->clear_image();
			ui.Detect_view->update();
			return;
		}
	}

	ui.Detect_view->cam_ratio = (double)cur_img.cols / (double)cur_img.rows;
	ui.Detect_view->set_image(QImage((const unsigned char*)cur_img.data, cur_img.cols,
		cur_img.rows, cur_img.step, QImage::Format_BGR888).copy());
	ui.Detect_view->update();

}

void SFM_module::add_scale_lable()
{
	begin_scale_box.push_back(new QSpinBox);
	_scale_label.push_back(new QLabel);
	_2scale_label.push_back(new QLabel);
	_3scale_label.push_back(new QLabel);
	_4scale_label.push_back(new QLabel);
	end_scale_box.push_back(new QSpinBox);
	value_scale_box.push_back(new QDoubleSpinBox);
	begin_scale_box[begin_scale_box.size() - 1]->setMinimum(-1);
	begin_scale_box[begin_scale_box.size() - 1]->setMaximum(1000);
	end_scale_box[end_scale_box.size() - 1]->setMinimum(-1);
	end_scale_box[end_scale_box.size() - 1]->setMaximum(1000);
	value_scale_box[end_scale_box.size() - 1]->setDecimals(6);
	value_scale_box[end_scale_box.size() - 1]->setMinimum(0.000001);
	value_scale_box[end_scale_box.size() - 1]->setMaximum(999999);
	begin_scale_box[begin_scale_box.size() - 1]->setValue(-1);
	end_scale_box[end_scale_box.size() - 1]->setValue(-1);
	value_scale_box[value_scale_box.size() - 1]->setValue(1);
	_scale_label[_scale_label.size() - 1]->setText(tr("-"));
	_2scale_label[_scale_label.size() - 1]->setText(tr(":"));
	_3scale_label[_scale_label.size() - 1]->setText(tr("mm"));
	_4scale_label[_scale_label.size() - 1]->setText(tr("Invalid"));
	QPalette pe;
	pe.setColor(QPalette::WindowText, Qt::darkRed);
	_4scale_label[_scale_label.size() - 1]->setPalette(pe);

	ui.gridLayout_scale->addWidget(begin_scale_box[begin_scale_box.size() - 1], begin_scale_box.size() - 1, 0, 1, 1);
	ui.gridLayout_scale->addWidget(_scale_label[_scale_label.size() - 1], _scale_label.size() - 1, 1, 1, 1);
	ui.gridLayout_scale->addWidget(end_scale_box[end_scale_box.size() - 1], end_scale_box.size() - 1, 2, 1, 1);
	ui.gridLayout_scale->addWidget(_2scale_label[_2scale_label.size() - 1], _2scale_label.size() - 1, 3, 1, 1);
	ui.gridLayout_scale->addWidget(value_scale_box[value_scale_box.size() - 1], value_scale_box.size() - 1, 4, 1, 1);
	ui.gridLayout_scale->addWidget(_3scale_label[_3scale_label.size() - 1], _3scale_label.size() - 1, 5, 1, 1);
	ui.gridLayout_scale->addWidget(_4scale_label[_4scale_label.size() - 1], _4scale_label.size() - 1, 6, 1, 1);

	connect(begin_scale_box[begin_scale_box.size() - 1], SIGNAL(valueChanged(int)), this, SLOT(test_valid_scale()));
	connect(end_scale_box[end_scale_box.size() - 1], SIGNAL(valueChanged(int)), this, SLOT(test_valid_scale()));
	connect(value_scale_box[value_scale_box.size() - 1], SIGNAL(valueChanged(double)), this, SLOT(test_valid_scale()));
}

void SFM_module::clear_scale_lable()
{
	for (unsigned short ii = 0; ii < begin_scale_box.size(); ii++)
	{
		ui.gridLayout_scale->removeWidget(begin_scale_box[ii]);
		ui.gridLayout_scale->removeWidget(_scale_label[ii]);
		ui.gridLayout_scale->removeWidget(_2scale_label[ii]);
		ui.gridLayout_scale->removeWidget(_3scale_label[ii]);
		ui.gridLayout_scale->removeWidget(_4scale_label[ii]);
		ui.gridLayout_scale->removeWidget(end_scale_box[ii]);
		ui.gridLayout_scale->removeWidget(value_scale_box[ii]);
		disconnect(begin_scale_box[ii]);
		disconnect(end_scale_box[ii]);
		disconnect(value_scale_box[ii]);
		delete begin_scale_box[ii];
		delete _scale_label[ii];
		delete _2scale_label[ii];
		delete _3scale_label[ii];
		delete _4scale_label[ii];
		delete end_scale_box[ii];
		delete value_scale_box[ii];
	}
	begin_scale_box.clear();
	_scale_label.clear();
	_2scale_label.clear();
	_3scale_label.clear();
	_4scale_label.clear();
	end_scale_box.clear();
	value_scale_box.clear();
	test_valid_scale();
}

void SFM_module::test_valid_scale()
{
	int vliad_number = 0;
	double sum_real_L = 0;
	double sum_cal_L = 0;
	for (int ii = 0; ii < begin_scale_box.size(); ii++)
	{
		if (begin_scale_box[ii]->value() == end_scale_box[ii]->value())
		{
			_4scale_label[ii]->setText(tr("Invalid"));
			QPalette pe;
			pe.setColor(QPalette::WindowText, Qt::darkRed);
			_4scale_label[ii]->setPalette(pe);
			continue;
		}
		int begin_index, end_index;
		bool has_found = false;
		for (int jj = 0; jj < Result_xyz_serial.size(); jj++)
		{
			if (Result_xyz_serial[jj].code == begin_scale_box[ii]->value())
			{
				begin_index = jj;
				has_found = true;
				break;
			}
		}
		if (!has_found)
		{
			_4scale_label[ii]->setText(tr("Invalid"));
			QPalette pe;
			pe.setColor(QPalette::WindowText, Qt::darkRed);
			_4scale_label[ii]->setPalette(pe);
			continue;
		}
		has_found = false;
		for (int jj = 0; jj < Result_xyz_serial.size(); jj++)
		{
			if (Result_xyz_serial[jj].code == end_scale_box[ii]->value())
			{
				end_index = jj;
				has_found = true;
				break;
			}
		}
		if (!has_found)
		{
			_4scale_label[ii]->setText(tr("Invalid"));
			QPalette pe;
			pe.setColor(QPalette::WindowText, Qt::darkRed);
			_4scale_label[ii]->setPalette(pe);
			continue;
		}

		_4scale_label[ii]->setText(tr("Pass"));
		QPalette pe;
		pe.setColor(QPalette::WindowText, Qt::darkGreen);
		_4scale_label[ii]->setPalette(pe);
		
		sum_cal_L += sqrt(pow(Result_xyz_serial[begin_index].x - Result_xyz_serial[end_index].x, 2) +
			pow(Result_xyz_serial[begin_index].y - Result_xyz_serial[end_index].y, 2) +
			pow(Result_xyz_serial[begin_index].z - Result_xyz_serial[end_index].z, 2));
		sum_real_L += value_scale_box[ii]->value();
		vliad_number++;
	}
	if (vliad_number != 0)
	{
		if (scale_for_sfm != sum_real_L / sum_cal_L)
		{
			scale_for_sfm = sum_real_L / sum_cal_L;
			update_points_view();
		}
	}
	else
	{
		if (scale_for_sfm != 1)
		{
			scale_for_sfm = 1;
			update_points_view();
		}
	}
}
void SFM_module::show_draw_setting_window()
{
	PCL_setting_window->show();
}

void SFM_module::show_setting_window()
{
	SFM_setting_window->show();
}
void SFM_module::setting_window_change()
{
		
	need_recal_keypoints = true;
}
void SFM_module::pcl_setting_window_change()
{
	update_points_view();
}

void SFM_module::open_images()
{
	//sbsbsb
	//QString  ImagePath = "C:\\Users\\Yzy_seu\\Desktop\\Image\\4";
	//QDir dir(ImagePath);
	//QStringList ImageList;
	//ImageList << "*.bmp";
	//dir.setNameFilters(ImageList);
	//SFM_setting_window->ui.r_k_doubleSpinBox->setValue(1);
	//SFM_setting_window->ui.r_k_1_doubleSpinBox->setValue(2);
	//SFM_setting_window->ui.r_k_2_doubleSpinBox->setValue(3);
	//SFM_setting_window->ui.BA_interval_spinBox->setValue(1);
	//SFM_setting_window->ui.min_radius_doubleSpinBox->setValue(10);
	//SFM_setting_window->ui.max_aspect_ratio_doubleSpinBox->setValue(4);
	//SFM_setting_window->ui.min_view_spinBox->setValue(2);
	//SFM_setting_window->ui.uniform_f_checkBox->setChecked(true);
	//SFM_setting_window->ui.optimize_checkBox->setChecked(true);
	//SFM_setting_window->ui.fixed_camera_checkBox->setChecked(true);
	//ui.reproject_error_limit_doubleSpinBox->setValue(0.25);
	////SFM_setting_window->ui.unlimit_checkBox->setChecked(false);
	////SFM_setting_window->ui.sub_pixel_method_comboBox->setCurrentIndex(0);
	//ui.focal_pixel_comboBox->setCurrentIndex(3);

	//QFile file("C:\\Users\\Yzy_seu\\Desktop\\Image\\read2.csv");
	//QStringList lines;/*行数据*/

	//std::vector<std::vector<int>> s_index;
	//if (file.open(QIODevice::ReadOnly))
	//{
	//	QTextStream stream_text(&file);
	//	while (!stream_text.atEnd())
	//	{
	//		lines.push_back(stream_text.readLine());
	//	}
	//	for (int j = 0; j < lines.size(); j++)
	//	{
	//		QString line = lines.at(j);
	//		QStringList split = line.split(",");/*列数据*/
	//		std::vector<int> s_now;
	//		for (int col = 0; col < split.size(); col++)
	//		{
	//			s_now.push_back(split[col].toInt());
	//		}
	//		s_index.push_back(s_now);
	//	}
	//	file.close();
	//}


	//QFile file2("C:\\Users\\Yzy_seu\\Desktop\\Image\\read4.csv");
	//QStringList lines2;/*行数据*/
	//R_save.clear();
	//T_save.clear();
	//if (file2.open(QIODevice::ReadOnly))
	//{
	//	QTextStream stream_text(&file2);
	//	while (!stream_text.atEnd())
	//	{
	//		lines2.push_back(stream_text.readLine());
	//	}
	//	for (int j = 0; j < lines2.size(); j++)
	//	{
	//		QString line = lines2.at(j);
	//		QStringList split = line.split(",");/*列数据*/
	//		std::vector<double> s_now;
	//		for (int col = 0; col < split.size(); col++)
	//		{
	//			s_now.push_back(split[col].toDouble());
	//		}
	//		cv::Mat R_temp = cv::Mat::ones(3, 3, CV_64F);
	//		cv::Mat T_temp = cv::Mat::ones(3, 1, CV_64F);
	//		T_temp.at<double>(1, 0) = s_now.at(0);
	//		T_temp.at<double>(0, 0) = s_now.at(1);
	//		T_temp.at<double>(2, 0) = s_now.at(2);
	//		R_temp.at<double>(1, 0) = s_now.at(3);
	//		R_temp.at<double>(0, 0) = s_now.at(4);
	//		R_temp.at<double>(2, 0) = s_now.at(5);
	//		R_temp.at<double>(1, 1) = s_now.at(6);
	//		R_temp.at<double>(0, 1) = s_now.at(7);
	//		R_temp.at<double>(2, 1) = s_now.at(8);
	//		R_temp.at<double>(1, 2) = -s_now.at(9);
	//		R_temp.at<double>(0, 2) = -s_now.at(10);
	//		R_temp.at<double>(2, 2) = -s_now.at(11);
	//		R_save.push_back(R_temp);
	//		T_save.push_back(T_temp);
	//	}
	//	file2.close();
	//}

	//SFM_setting_window->ui.circle_correctin_checkBox->setChecked(true);
	//SFM_setting_window->ui.circle_iter_spinBox->setValue(3);
	//for (int sample = 2; sample <=2; sample++)
	//{
	//	for (int iter = 0; iter < s_index.size(); iter++)
	//	{
	//		std::cout << "sample=" << sample << ";iter=" << iter << std::endl; 
	//		QFileInfo fi("C:\\Users\\Yzy_seu\\Desktop\\Image\\Dual_after\\Result_" + QString("%1").arg(sample, 2, 10, QLatin1Char('0')) + "_" + QString("%1").arg(iter, 4, 10, QLatin1Char('0')) + ".csv");
	//		if (fi.exists())
	//		{
	//			continue;
	//		}
	//		clear_all_data();
	//		std::vector<QString> dir_sam;
	//		for (int ss = 0; ss < sample;ss++)
	//		{
	//			dir_sam.push_back("C:\\Users\\Yzy_seu\\Desktop\\Image\\4\\" + dir[s_index[iter][ss]]);

	//		}
	//		index_for_L = s_index[iter][0];
	//		index_for_R = s_index[iter][1];
	//		for (int ii = 0; ii < dir_sam.size(); ii++)
	//		{
	//			Image_serial_name.push_back(dir_sam[ii]);
	//			KeyPoint_Enable_serial.push_back(true);
	//		}
	//		QDir file_temp = QDir(dir_sam[0]);
	//		file_temp.cdUp();
	//		ui.open_path_lineEdit->setText(file_temp.absolutePath());
	//		update_table_view();
	//		if (Image_serial_name.size() == 0)
	//		{
	//			ui.show_detect_spinBox->setEnabled(false);
	//			ui.show_detect_horizontalSlider->setEnabled(false);
	//			ui.show_detect_spinBox->setMinimum(1);
	//			ui.show_detect_spinBox->setMaximum(1);
	//			ui.show_detect_horizontalSlider->setMinimum(1);
	//			ui.show_detect_horizontalSlider->setMaximum(1);
	//			ui.show_detect_spinBox->setValue(1);
	//			ui.show_detect_horizontalSlider->setValue(1);
	//			ui.Detect_view->clear_image();
	//		}
	//		else
	//		{
	//			ui.show_detect_spinBox->setEnabled(true);
	//			ui.show_detect_horizontalSlider->setEnabled(true);
	//			ui.show_detect_spinBox->setMinimum(1);
	//			ui.show_detect_spinBox->setMaximum(Image_serial_name.size());
	//			ui.show_detect_horizontalSlider->setMinimum(1);
	//			ui.show_detect_horizontalSlider->setMaximum(Image_serial_name.size());
	//			ui.show_detect_spinBox->setValue(1);
	//			ui.show_detect_horizontalSlider->setValue(1);
	//		}
	//		update_show_view();
	//		need_recal_keypoints = true;
	//		calculate_multiview();
	// 
	//		if (Result_xyz_serial.size() != 0 && Result_err_every.size()!=0)
	//		{
	//			update_points_view();
	//			save_SFMResult_file_csv("C:\\Users\\Yzy_seu\\Desktop\\Image\\Dual_after\\Result_" + QString("%1").arg(sample, 2, 10, QLatin1Char('0')) + "_" + QString("%1").arg(iter, 4, 10, QLatin1Char('0')) + ".csv");
	//		}
	//	}
	//}
	if (ui.data_image_radioButton->isChecked())
	{
		QString strs;
		QStringList file_list, output_name;
		QStringList str_path_list = QFileDialog::getOpenFileNames(this, tr("Select Images"), tr(""), tr("Image(*.bmp *.jpg *.png *.tiff *.jfif);;All(*)"));
		if (str_path_list.size() == 0)
		{
			return;
		}
		clear_all_data();


		for (int ii = 0; ii < str_path_list.size(); ii++)
		{
			Image_serial_name.push_back(str_path_list[ii]);
			KeyPoint_Enable_serial.push_back(true);
		}
		QDir file_temp = QDir(str_path_list[0]);
		file_temp.cdUp();
		ui.open_path_lineEdit->setText(file_temp.absolutePath());
		update_table_view();
		if (Image_serial_name.size() == 0)
		{
			ui.show_detect_spinBox->setEnabled(false);
			ui.show_detect_horizontalSlider->setEnabled(false);
			ui.show_detect_spinBox->setMinimum(1);
			ui.show_detect_spinBox->setMaximum(1);
			ui.show_detect_horizontalSlider->setMinimum(1);
			ui.show_detect_horizontalSlider->setMaximum(1);
			ui.show_detect_spinBox->setValue(1);
			ui.show_detect_horizontalSlider->setValue(1);
			ui.Detect_view->clear_image();
		}
		else
		{
			ui.show_detect_spinBox->setEnabled(true);
			ui.show_detect_horizontalSlider->setEnabled(true);
			ui.show_detect_spinBox->setMinimum(1);
			ui.show_detect_spinBox->setMaximum(Image_serial_name.size());
			ui.show_detect_horizontalSlider->setMinimum(1);
			ui.show_detect_horizontalSlider->setMaximum(Image_serial_name.size());
			ui.show_detect_spinBox->setValue(1);
			ui.show_detect_horizontalSlider->setValue(1);
		}
		update_show_view();
		need_recal_keypoints = true;
	}
	else
	{
		select_path_read_keypoint();
	}
}

void SFM_module::update_camera_params()
{
	QSettings settings;
	QStringList SFM_camera = settings.value("SFM_camera").toStringList();
	QStringList SFM_camera_focal = settings.value("SFM_camera_focal").toStringList();
	QStringList SFM_camera_pixel = settings.value("SFM_camera_pixel").toStringList();
	ui.focal_pixel_comboBox->clear();
	SFM_camera_focal_vector.clear();
	SFM_camera_pixel_vector.clear();
	SFM_camera_name_vector.clear();
	ui.focal_pixel_comboBox->addItem(tr("None"));
	ui.focal_pixel_comboBox->addItem(tr("User"));
	for (int ii = 0; ii < SFM_camera.size(); ii++)
	{
		if (ii >= SFM_camera_focal.size() || ii >= SFM_camera_pixel.size())
		{
			continue;
		}
		ui.focal_pixel_comboBox->addItem(SFM_camera[ii]);
		SFM_camera_name_vector.push_back(SFM_camera[ii]);
		SFM_camera_focal_vector.push_back(SFM_camera_focal[ii].toDouble());
		SFM_camera_pixel_vector.push_back(SFM_camera_pixel[ii].toDouble());
	}
}
void SFM_module::update_camera_params_enable()
{
	if (ui.focal_pixel_comboBox->count() <= 0)
	{
		return;
	}
	ui.save_focal_pixel_pushButton->setText(tr("Save"));
	if (ui.focal_pixel_comboBox->currentIndex() == 0)
	{
		ui.focal_length_doubleSpinBox->setEnabled(false);
		ui.pixel_size_doubleSpinBox->setEnabled(false);
		ui.save_focal_pixel_pushButton->setEnabled(false);
		ui.delete_focal_pixel_pushButton->setEnabled(false);
	}
	else if (ui.focal_pixel_comboBox->currentIndex() == 1)
	{
		ui.focal_length_doubleSpinBox->setEnabled(true);
		ui.pixel_size_doubleSpinBox->setEnabled(true);
		ui.save_focal_pixel_pushButton->setEnabled(true);
		ui.delete_focal_pixel_pushButton->setEnabled(false);
	}
	else
	{
		ui.focal_length_doubleSpinBox->setValue(SFM_camera_focal_vector[ui.focal_pixel_comboBox->currentIndex() - 2]);
		ui.pixel_size_doubleSpinBox->setValue(SFM_camera_pixel_vector[ui.focal_pixel_comboBox->currentIndex() - 2]);
		ui.focal_length_doubleSpinBox->setEnabled(true);
		ui.pixel_size_doubleSpinBox->setEnabled(true);
		ui.save_focal_pixel_pushButton->setEnabled(true);
		ui.delete_focal_pixel_pushButton->setEnabled(true);
		ui.save_focal_pixel_pushButton->setText(tr("Update"));
	}
}
void SFM_module::save_camera_params()
{
	if (ui.focal_pixel_comboBox->currentIndex() == 1)
	{
		bool bOk = false;
		QString sName = QInputDialog::getText(this,
			tr("SetCameraFlag"),
			tr("Camera Name:"),
			QLineEdit::Normal,
			tr(""),
			&bOk
		);
		if (!bOk)
		{
			return;
		}
		sName = sName + tr("-F") + QString::number(ui.focal_length_doubleSpinBox->value(), 'f', 2) 
			+ tr("-S") + QString::number(ui.pixel_size_doubleSpinBox->value(), 'f', 2);
		bool is_exist = false;
		int cur_index;
		for (int ii = 0; ii < SFM_camera_name_vector.size(); ii++)
		{
			if (!SFM_camera_name_vector[ii].compare(sName))
			{
				is_exist = true;
				cur_index = ii;
				break;
			}
		}
		if (is_exist)
		{
			if ((cur_index + 2) < ui.focal_pixel_comboBox->count())
			{
				ui.focal_pixel_comboBox->setCurrentIndex(cur_index + 2);
			}
			QMessageBox::warning(this, tr("Save camera parameters"), tr("Camera already exists!"), QMessageBox::Ok, QMessageBox::Ok);
			return;
		}
		SFM_camera_name_vector.push_back(sName);
		SFM_camera_focal_vector.push_back(ui.focal_length_doubleSpinBox->value());
		SFM_camera_pixel_vector.push_back(ui.pixel_size_doubleSpinBox->value());

		QStringList SFM_camera;
		QStringList SFM_camera_focal;
		QStringList SFM_camera_pixel;
		for (int ii = 0; ii < SFM_camera_name_vector.size(); ii++)
		{
			SFM_camera.push_back(SFM_camera_name_vector[ii]);
			SFM_camera_focal.push_back(QString::number(SFM_camera_focal_vector[ii], 'f', 10));
			SFM_camera_pixel.push_back(QString::number(SFM_camera_pixel_vector[ii], 'f', 10));
		}

		QSettings settings;
		settings.setValue("SFM_camera", SFM_camera);
		settings.setValue("SFM_camera_focal", SFM_camera_focal);
		settings.setValue("SFM_camera_pixel", SFM_camera_pixel);
		settings.sync();
		update_camera_params();
		ui.focal_pixel_comboBox->setCurrentIndex(SFM_camera.size() + 1);
	}
	else if (ui.focal_pixel_comboBox->currentIndex() > 1)
	{
		int index = ui.focal_pixel_comboBox->currentIndex() - 2;
		if (index < 0 || index >= SFM_camera_name_vector.size())
		{
			return;
		}
		SFM_camera_focal_vector[index] = ui.focal_length_doubleSpinBox->value();
		SFM_camera_pixel_vector[index] = ui.pixel_size_doubleSpinBox->value();

		QStringList SFM_camera;
		QStringList SFM_camera_focal;
		QStringList SFM_camera_pixel;
		for (int ii = 0; ii < SFM_camera_name_vector.size(); ii++)
		{
			SFM_camera.push_back(SFM_camera_name_vector[ii]);
			SFM_camera_focal.push_back(QString::number(SFM_camera_focal_vector[ii], 'f', 10));
			SFM_camera_pixel.push_back(QString::number(SFM_camera_pixel_vector[ii], 'f', 10));
		}

		QSettings settings;
		settings.setValue("SFM_camera", SFM_camera);
		settings.setValue("SFM_camera_focal", SFM_camera_focal);
		settings.setValue("SFM_camera_pixel", SFM_camera_pixel);
		settings.sync();
		update_camera_params();
		ui.focal_pixel_comboBox->setCurrentIndex(SFM_camera.size() + 1);
	}
}
void SFM_module::delete_camera_params()
{
	int index = ui.focal_pixel_comboBox->currentIndex() - 2;
	if (index < 0 || index >= SFM_camera_name_vector.size())
	{
		return;
	}
	QMessageBox::StandardButton box;
	box = QMessageBox::question(this, tr("Confirm"), tr("Delete?"), QMessageBox::Yes | QMessageBox::No);
	if (box == QMessageBox::No)
	{
		return;
	}

	SFM_camera_name_vector.erase(SFM_camera_name_vector.begin() + index);
	SFM_camera_focal_vector.erase(SFM_camera_focal_vector.begin() + index);
	SFM_camera_pixel_vector.erase(SFM_camera_pixel_vector.begin() + index);

	QStringList SFM_camera;
	QStringList SFM_camera_focal;
	QStringList SFM_camera_pixel;
	for (int ii = 0; ii < SFM_camera_name_vector.size(); ii++)
	{
		SFM_camera.push_back(SFM_camera_name_vector[ii]);
		SFM_camera_focal.push_back(QString::number(SFM_camera_focal_vector[ii], 'f', 10));
		SFM_camera_pixel.push_back(QString::number(SFM_camera_pixel_vector[ii], 'f', 10));
	}

	QSettings settings;
	settings.setValue("SFM_camera", SFM_camera);
	settings.setValue("SFM_camera_focal", SFM_camera_focal);
	settings.setValue("SFM_camera_pixel", SFM_camera_pixel);
	settings.sync();
	update_camera_params();
}
void SFM_module::clear_all_data()
{
	Image_serial_name.clear();
	Result_err_every.clear();
	Result_xyz_serial.clear();
	Result_xyz_after_op_serial.clear();
	Result_R_after_op.clear();
	Result_T_after_op.clear();
	Result_Ori_after_op.clear();
	xyz_need_levelling_serial.clear();
	xyz_need_levelling_code_serial.clear();
	KeyPoint_code_serial.clear();
	KeyPoint_code_serial_update.clear();
	ImageWidth_serial.clear();
	ImageHeight_serial.clear();
	KeyPoint_Enable_serial.clear();
	Result_contours_serial.clear();
	Result_R.clear();
	Result_T.clear();
	Result_Ori.clear();
	Keypoints_in_queue.clear();
	Keypoints_contours_in_queue.clear();
	Keypoints_index_in_queue.clear();
	Keypoints_R_in_queue.clear();
	Keypoints_T_in_queue.clear();
	update_table_view();
	update_show_view();
	clear_scale_lable();
	ui.result_fx_doubleSpinBox->setValue(0);
	ui.result_fy_doubleSpinBox->setValue(0);
	ui.result_cx_doubleSpinBox->setValue(0);
	ui.result_cy_doubleSpinBox->setValue(0);
	ui.result_shear_doubleSpinBox->setValue(0);
	ui.result_dis_K_doubleSpinBox_1->setValue(0);
	ui.result_dis_K_doubleSpinBox_2->setValue(0);
	ui.result_dis_K_doubleSpinBox_3->setValue(0);
	ui.result_dis_K_doubleSpinBox_4->setValue(0);
	ui.result_dis_K_doubleSpinBox_5->setValue(0);
	ui.result_dis_K_doubleSpinBox_6->setValue(0);
	ui.result_dis_P_doubleSpinBox->setValue(0);
	ui.result_dis_P_doubleSpinBox_2->setValue(0);
	ui.result_dis_T_doubleSpinBox_1->setValue(0);
	ui.result_dis_T_doubleSpinBox_2->setValue(0);
	ui.result_dis_T_doubleSpinBox_3->setValue(0);
	ui.result_dis_T_doubleSpinBox_4->setValue(0);

	need_recal_keypoints = true;
}

void SFM_module::loop_sleep(int msec)
{
	QTime dieTime = QTime::currentTime().addMSecs(msec);

	while (QTime::currentTime() < dieTime) {
		QCoreApplication::processEvents(QEventLoop::AllEvents, 100);
	}
}

void SFM_module::calculate_multiview()
{

	if (!ui.focal_length_doubleSpinBox->isEnabled())
	{
		QMessageBox::warning(this, tr("Parameter check"), tr("No focal length and pixel information, please set!"), QMessageBox::Ok, QMessageBox::Ok);
		return;
	}
	cam_log_ui->show();
	cam_log_ui->setWindowOpacity(0.5);
	cam_log_ui->ui.textEdit->clear();
	QString str_inf_for_log = "";
	if (Image_serial_name.size() == 0)
	{
		QMessageBox::warning(this, tr("Detect Circles"), tr("No images!"), QMessageBox::Ok, QMessageBox::Ok);
		return;
	}
	ui.view_tabWidget->setCurrentIndex(0);
	Result_xyz_serial.clear();
	Result_xyz_after_op_serial.clear();
	Result_R_after_op.clear();
	Result_T_after_op.clear();
	Result_Ori_after_op.clear();
	xyz_need_levelling_serial.clear();
	xyz_need_levelling_code_serial.clear();
	Result_R.clear();
	Result_T.clear();
	Result_Ori.clear();
	Result_err_every.clear();
	Result_err_every.clear();
	DetectContoursMethod contour_method;
	switch (SFM_setting_window->ui.feature_detect_method_comboBox->currentIndex())
	{
	case 0:
		contour_method = DetectContoursMethod::Sobel_Method;
		break;
	case 1:
		contour_method = DetectContoursMethod::CANNY_Method;
		break;
	case 2:
		contour_method = DetectContoursMethod::Prewitt_Method;
		break;
	case 3:
		contour_method = DetectContoursMethod::Roberts_Method;
		break;
	case 4:
		contour_method = DetectContoursMethod::OTSU_Method;
		break;
	case 5:
		contour_method = DetectContoursMethod::ADAPTIVE_THRESH_Method;
		break;
	case 6:
		contour_method = DetectContoursMethod::Block_OTSU;
		break;
	default:
		contour_method = DetectContoursMethod::Sobel_Method;
		break;
	}
	SubPixelPosMethod sub_pixel_method;
	switch (SFM_setting_window->ui.sub_pixel_method_comboBox->currentIndex())
	{
	case 0:
		sub_pixel_method = SubPixelPosMethod::NoSubPixel_Match;
		break;
	case 1:
		sub_pixel_method = SubPixelPosMethod::Interpolation_Ellipse_Match;
		break;
	case 2:
		sub_pixel_method = SubPixelPosMethod::Interpolation_Rotaion_Ellipse_Match;
		break;
	case 3:
		sub_pixel_method = SubPixelPosMethod::Surface_fit_Ellipse_Match;
		break;
	case 4:
		sub_pixel_method = SubPixelPosMethod::Gauss_surface_fit_Ellipse_Match;
		break;
	case 5:
		sub_pixel_method = SubPixelPosMethod::Gauss_Curve_Fit;
		break;
	case 6:
		sub_pixel_method = SubPixelPosMethod::Gray_Moment;
		break;
	case 7:
		sub_pixel_method = SubPixelPosMethod::Binary_Centroid;
		break;
	case 8:
		sub_pixel_method = SubPixelPosMethod::Gray_Centroid;
		break;
	case 9:
		sub_pixel_method = SubPixelPosMethod::Squared_Gray_Centroid;
		break;
	default:
		sub_pixel_method = SubPixelPosMethod::Gray_Centroid;
		break;
	}

	MarkPointColorType mark_type;
	switch (SFM_setting_window->ui.color_comboBox->currentIndex())
	{
	case 0:
		mark_type = MarkPointColorType::BlackDownWhiteUp;
		break;
	case 1:
		mark_type = MarkPointColorType::WhiteDownBlackUp;
		break;
	case 2:
		mark_type = MarkPointColorType::Uncertainty;
		break;
	default:
		mark_type = MarkPointColorType::BlackDownWhiteUp;
		break;
	}
	CodePointBitesType circle_code_type;
	switch (SFM_setting_window->ui.code_bits_comboBox->currentIndex())
	{
	case 0:
		circle_code_type = CodePointBitesType::CodeBites15;
		break;
	case 1:
		circle_code_type = CodePointBitesType::CodeBites12;
		break;
	default:
		circle_code_type = CodePointBitesType::CodeBites15;
		break;
	}
	if (need_recal_keypoints)
	{
		str_inf_for_log = tr("----------------------------------------") + "\n" + 
			tr("Need to identify the encoding circle.") + "\n" + str_inf_for_log;
		cam_log_ui->ui.textEdit->setText(str_inf_for_log);
		ui.Detect_view->OnPresetImage_from_software();
		ui.show_detect_horizontalSlider->setValue(1);
		ui.show_detect_spinBox->setValue(1);

		in_calculation = true;
		KeyPoint_code_serial.clear();
		KeyPoint_code_serial_update.clear();
		ImageWidth_serial.clear();
		ImageHeight_serial.clear();
		Result_contours_serial.clear();

		QThreadPool::globalInstance()->setMaxThreadCount(SFM_setting_window->ui.thread_spinBox->value());



		std::vector<int> thread_for_detect_calcu_index;
		std::vector<bool> thread_for_detect_calcu_finished;
		std::vector<deltect_coded_circle_thread*> thread_for_detect_calcu;
		for (int ii = 0; ii < Image_serial_name.size(); ii++)
		{
			std::vector<std::vector<cv::Point>> contours_now;
			std::vector<Coded_detect_inf> dect_code_temp;
			KeyPoint_code_serial.push_back(dect_code_temp);
			Result_contours_serial.push_back(contours_now);
			ImageWidth_serial.push_back(0);
			ImageHeight_serial.push_back(0);
			if (!KeyPoint_Enable_serial[ii])
			{
				continue;
			}
			thread_for_detect_calcu_index.push_back(ii);
			thread_for_detect_calcu_finished.push_back(false);
			thread_for_detect_calcu.push_back(new deltect_coded_circle_thread(Image_serial_name[ii], SFM_setting_window->ui.r_k_doubleSpinBox->value(),
				SFM_setting_window->ui.r_k_1_doubleSpinBox->value(), SFM_setting_window->ui.r_k_2_doubleSpinBox->value(),
				SFM_setting_window->ui.min_radius_doubleSpinBox->value(), SFM_setting_window->ui.max_radius_doubleSpinBox->value(),
				SFM_setting_window->ui.max_radius_error_doubleSpinBox->value(), mark_type, circle_code_type,
				contour_method, sub_pixel_method, SFM_setting_window->ui.max_aspect_ratio_doubleSpinBox->value()
				, SFM_setting_window->ui.min_contour_points_spinBox->value(), SFM_setting_window->ui.min_contour_number_spinBox->value()
				, SFM_setting_window->ui.min_gap_spinBox->value()
				, SFM_setting_window->ui.max_foreground_std_spinBox->value(), SFM_setting_window->ui.max_background_std_spinBox->value()));
			thread_for_detect_calcu[thread_for_detect_calcu.size() - 1]->setAutoDelete(false);
		}

		QProgressDialog oQProgressDialog;
		oQProgressDialog.setWindowModality(Qt::ApplicationModal);
		oQProgressDialog.setWindowTitle(tr("Detect coded circle points..."));
		oQProgressDialog.setRange(0, thread_for_detect_calcu.size() + 1);
		oQProgressDialog.setMinimumDuration(0);

		for (int ii = 0; ii < thread_for_detect_calcu.size(); ii++)
		{
			QThreadPool::globalInstance()->start(thread_for_detect_calcu[ii]);
		}

		oQProgressDialog.setValue(1);
		qApp->processEvents();

		while (1)
		{
			loop_sleep(10);
			int finished_number = 0;
			if (oQProgressDialog.wasCanceled())
			{
				KeyPoint_code_serial.clear();
				KeyPoint_code_serial_update.clear();
				ImageWidth_serial.clear();
				ImageHeight_serial.clear();
				Result_contours_serial.clear();
				QMessageBox::information(this, tr("Process"), tr("User Stop"), QMessageBox::Ok);
				in_calculation = false;
				return;
			}
			for (int ii = 0; ii < thread_for_detect_calcu.size(); ii++)
			{
				if (thread_for_detect_calcu_finished[ii])
				{
					finished_number++;
					continue;
				}
				if (thread_for_detect_calcu[ii]->is_finish)
				{
					finished_number++;
					thread_for_detect_calcu_finished[ii] = true;
					cv::Mat code_point_temp = thread_for_detect_calcu[ii]->code_point_mat;
					if (thread_for_detect_calcu[ii]->success_cal && code_point_temp.rows > 0)
					{
						str_inf_for_log = tr("Image-") + QString::number(ii) +
							tr(" recognition successful, number of encoding points: ") + QString::number(code_point_temp.rows) +
							tr(".") + "\n" + str_inf_for_log;
						cam_log_ui->ui.textEdit->setText(str_inf_for_log);
						std::vector<Coded_detect_inf> dect_code_temp;
						for (unsigned int jj = 0; jj < code_point_temp.rows; jj++)
						{
							Coded_detect_inf TEMP_INF;
							TEMP_INF.code_num = code_point_temp.at<float>(jj, 0);
							TEMP_INF.x = code_point_temp.at<float>(jj, 1);
							TEMP_INF.y = code_point_temp.at<float>(jj, 2);
							TEMP_INF.fit_err = code_point_temp.at<float>(jj, 3);
							TEMP_INF.r_a = code_point_temp.at<float>(jj, 4);
							TEMP_INF.r_b = code_point_temp.at<float>(jj, 5);
							TEMP_INF.angle_in_pi = code_point_temp.at<float>(jj, 6);
							TEMP_INF.image_index = thread_for_detect_calcu_index[ii];
							dect_code_temp.push_back(TEMP_INF);
						}
						KeyPoint_code_serial[thread_for_detect_calcu_index[ii]] = dect_code_temp;
						Result_contours_serial[thread_for_detect_calcu_index[ii]] = thread_for_detect_calcu[ii]->contours_pixel;
						ImageWidth_serial[thread_for_detect_calcu_index[ii]] = (thread_for_detect_calcu[ii]->ori_image_thread.cols);
						ImageHeight_serial[thread_for_detect_calcu_index[ii]] = (thread_for_detect_calcu[ii]->ori_image_thread.rows);
						if (thread_for_detect_calcu_index[ii] + 1 >= ui.show_detect_horizontalSlider->value())
						{
							ui.show_detect_horizontalSlider->setValue(thread_for_detect_calcu_index[ii] + 1);
							ui.show_detect_spinBox->setValue(thread_for_detect_calcu_index[ii] + 1);
							update_show_view_once(thread_for_detect_calcu[ii]->ori_image_thread, thread_for_detect_calcu_index[ii]);
						}
					}
					else
					{
						str_inf_for_log = tr("Image-") + QString::number(ii) +
							tr(" recognition failed. ") + "\n" + str_inf_for_log;
						cam_log_ui->ui.textEdit->setText(str_inf_for_log);
						//KeyPoint_Enable_serial[thread_for_detect_calcu_index[ii]] = false;
					}
					delete thread_for_detect_calcu[ii];
				}
			}

			oQProgressDialog.setLabelText(QString::number(finished_number) + tr(" images have been detected."));
			oQProgressDialog.setValue(finished_number + 1);
			qApp->processEvents();
			if (finished_number == thread_for_detect_calcu.size())
			{
				break;
			}
		}

		update_show_view();
		need_recal_keypoints = false;
		oQProgressDialog.cancel();
		in_calculation = false;
	}

	KeyPoint_code_serial_update = KeyPoint_code_serial;
	if (KeyPoint_code_serial.size() < 2)
	{
		QMessageBox::warning(this, tr("SFM"), tr("The matching images less than 2!"), QMessageBox::Ok, QMessageBox::Ok);
		return;
	}


	str_inf_for_log = tr("----------------------------------------") + "\n" +
		tr("Pre-orientation:Initial Frame Pair Selection.") + "\n" + str_inf_for_log;
	cam_log_ui->ui.textEdit->setText(str_inf_for_log);
	qApp->processEvents();

	std::vector<std::vector<Coded_detect_inf>> KeyPoint_code_serial_copy;
	std::vector<std::vector<std::vector<cv::Point>> > KeyPoint_contour_serial_copy;
	std::vector<int> KeyPoint_code_serial_index_copy;
	for (int ii = 0; ii < KeyPoint_Enable_serial.size(); ii++)
	{
		if (KeyPoint_Enable_serial[ii])
		{
			KeyPoint_code_serial_copy.push_back(KeyPoint_code_serial[ii]);
			KeyPoint_code_serial_index_copy.push_back(ii);
			KeyPoint_contour_serial_copy.push_back(Result_contours_serial[ii]);
		}
	}
	std::priority_queue<match_pair_group> match_queue;
	std::vector<match_pair_group> match_need_cal;
	int max_number = 0;
	for (int ii = 0; ii < KeyPoint_code_serial_copy.size(); ii++)
	{
		for (int jj = ii + 1; jj < KeyPoint_code_serial_copy.size(); jj++)
		{
			match_need_cal.push_back(match_pair_group(find_same_coded_point(KeyPoint_code_serial_copy[ii], KeyPoint_code_serial_copy[jj]), ii, jj));
			match_queue.push(match_pair_group(find_same_coded_point(KeyPoint_code_serial_copy[ii], KeyPoint_code_serial_copy[jj]), ii, jj));
		}
	}
	if (match_queue.top().matching_number < 5)
	{
		QMessageBox::warning(this, tr("SFM"), tr("The maximum number of matching points for any image pair are less than 5!"), QMessageBox::Ok, QMessageBox::Ok);
		return;
	}


	QProgressDialog oQProgressDialog;
	oQProgressDialog.setWindowModality(Qt::ApplicationModal);
	oQProgressDialog.setWindowTitle(tr("Calculation"));
	oQProgressDialog.setLabelText(tr("Incremental reconstruction..."));
	oQProgressDialog.setRange(0, KeyPoint_code_serial_index_copy.size() + 2);
	oQProgressDialog.setValue(0);
	oQProgressDialog.setMinimumDuration(0);
	double fx = 1, fy = 1, cx = ImageWidth_serial[match_queue.top().position_x] * 0.5, 
		cy = ImageHeight_serial[match_queue.top().position_y] * 0.5;
	if (ui.camera_param_hardware_radioButton->isChecked())
	{
		fx = ui.focal_length_doubleSpinBox->value() / ui.pixel_size_doubleSpinBox->value() * 1000.0;
		fy = fx;
	}
	else if (ui.camera_param_calibrated_radioButton->isChecked())
	{
		if (has_loaded_camera_params)
		{
			fx = loaded_camera_params[0];
			fy = loaded_camera_params[1];
			cx = loaded_camera_params[2];
			cy = loaded_camera_params[3];
		}
		else
		{
			oQProgressDialog.cancel();
			QMessageBox::warning(this, tr("SFM"), tr("No calibration result loaded!"), QMessageBox::Ok, QMessageBox::Ok);
			return;
		}
	}

	cv::Mat K_camera(cv::Matx33d(
		fx, 0, cx,
		0, fy, cy,
		0, 0, 1));
	//计算初状态

	//sbsbsb
	/*bool need_optimize_C;
	bool need_optimize_F;
	double* range_f = new double[3];
	double* range_c = new double[3];
	if (SFM_setting_window->ui.unlimit_checkBox->isChecked())
	{
		delete[] range_f;
		delete[] range_c;
		range_f = nullptr;
		range_c = nullptr;
	}
	else
	{
		range_f[0] = fx;
		range_f[1] = fy;
		range_f[2] = ((double)SFM_setting_window->ui.limit_spinBox->value()) / 100.0 * fx < ((double)SFM_setting_window->ui.limit_spinBox->value()) / 100.0 * fy
			? ((double)SFM_setting_window->ui.limit_spinBox->value()) / 100.0 * fx : ((double)SFM_setting_window->ui.limit_spinBox->value()) / 100.0 * fy;
		range_c[0] = cx;
		range_c[1] = cy;
		range_c[2] = ((double)SFM_setting_window->ui.limit_spinBox->value()) / 100.0 * cx < ((double)SFM_setting_window->ui.limit_spinBox->value()) / 100.0 * cy
			? ((double)SFM_setting_window->ui.limit_spinBox->value()) / 100.0 * cx : ((double)SFM_setting_window->ui.limit_spinBox->value()) / 100.0 * cy;
	}
	std::vector<sfm_3d_group> xyz_sfm_ori;
	double* K = new double[5];
	double* Dis_K = new double[6];
	double* Dis_P = new double[2];
	double* Dis_T = new double[4];
	K[0] = fx;
	K[1] = fy;
	K[2] = cx;
	K[3] = cy;
	K[4] = 0;
	Dis_K[0] = 0;
	Dis_K[1] = 0;
	Dis_K[2] = 0;
	Dis_K[3] = 0;
	Dis_K[4] = 0;
	Dis_K[5] = 0;
	Dis_P[0] = 0;
	Dis_P[1] = 0;
	Dis_T[0] = 0;
	Dis_T[1] = 0;
	Dis_T[2] = 0;
	Dis_T[3] = 0;
	cv::Mat dis_coff = cv::Mat::eye(1, 12, CV_64FC1);
	dis_coff.at<double>(0) = Dis_K[0];
	dis_coff.at<double>(1) = Dis_K[1];
	dis_coff.at<double>(2) = Dis_P[0];
	dis_coff.at<double>(3) = Dis_P[1];
	dis_coff.at<double>(4) = Dis_K[2];
	dis_coff.at<double>(5) = Dis_K[3];
	dis_coff.at<double>(6) = Dis_K[4];
	dis_coff.at<double>(7) = Dis_K[5];
	dis_coff.at<double>(8) = Dis_T[0];
	dis_coff.at<double>(9) = Dis_T[2];
	dis_coff.at<double>(10) = Dis_T[1];
	dis_coff.at<double>(11) = Dis_T[3];*/



	//sbsbsb
	std::vector<cv::Point2d> P1_ori_first, P2_ori_first;
	std::vector<int> Code_ori_first;
	std::vector<sfm_3d_group> xyz_sfm_ori;
	cv::Mat R_ori_first, T_ori_first;	//旋转矩阵和平移向量
	cv::Mat mask_ori_first;	//mask中大于零的点代表匹配点，等于零代表失配点
	while (match_queue.size())
	{
		if (match_queue.top().matching_number < 4)
		{
			break;
		}
		get_matched_points(KeyPoint_code_serial_copy[match_queue.top().position_x],
			KeyPoint_code_serial_copy[match_queue.top().position_y], P1_ori_first, P2_ori_first, Code_ori_first);
		if (cal_transform(K_camera, P1_ori_first, P2_ori_first, R_ori_first, T_ori_first, mask_ori_first, SFM_setting_window->ui.object_comboBox->currentIndex()) >= 4)
		{
			Eigen::Matrix<double, 3, 3> matrix_temp;
			cv::cv2eigen(R_ori_first, matrix_temp);
			Eigen::Quaterniond qua_temp = Eigen::Quaterniond(matrix_temp);
			qua_temp.normalize();
			double max_q_value = abs(qua_temp.x());
			max_q_value = max_q_value > abs(qua_temp.y()) ? max_q_value : abs(qua_temp.y());
			max_q_value = max_q_value > abs(qua_temp.z()) ? max_q_value : abs(qua_temp.z());
			max_q_value = max_q_value > abs(qua_temp.w()) ? max_q_value : abs(qua_temp.w());
			if (((abs(T_ori_first.at<double>(2)) / sqrt(T_ori_first.at<double>(0) * T_ori_first.at<double>(0) + T_ori_first.at<double>(1) * T_ori_first.at<double>(1) + T_ori_first.at<double>(2) * T_ori_first.at<double>(2)) < 0.8)
				|| max_q_value / sqrt(qua_temp.x() * qua_temp.x() + qua_temp.y() * qua_temp.y() + qua_temp.z() * qua_temp.z() + qua_temp.w() * qua_temp.w()) < 0.99))
			{
				break;
			}
		}
		match_queue.pop();
	}
	if (match_queue.size() == 0 || match_queue.top().matching_number < 4)
	{
		oQProgressDialog.cancel();
		QMessageBox::warning(this, tr("SFM"), tr("Failed to initialize the posture in the predetermined direction: First KeyFrame instability!"), QMessageBox::Ok, QMessageBox::Ok);
		cam_log_ui->hide();
		return;
	}
	oQProgressDialog.setValue(2);
	str_inf_for_log = tr("----------------------------------------") + "\n" +
		tr("Successfully selected initial frame pair: Image") + QString::number(KeyPoint_code_serial_index_copy[match_queue.top().position_x] + 1)
		+ tr("-Image") + QString::number(KeyPoint_code_serial_index_copy[match_queue.top().position_y] + 1) + ".\n" + str_inf_for_log;
	cam_log_ui->ui.textEdit->setText(str_inf_for_log);
	qApp->processEvents();
	maskout_points(P1_ori_first, P2_ori_first, Code_ori_first, mask_ori_first);
	if (P1_ori_first.size() < 4)
	{
		oQProgressDialog.cancel();
		QMessageBox::warning(this, tr("SFM"), tr("Failed to initialize the posture in the predetermined direction: First KeyFrame points less than 4!"), QMessageBox::Ok, QMessageBox::Ok);
		cam_log_ui->hide();
		return;
	}
	std::vector<sfm_3d_group> xyz_ori_first = reconstruct_3d_first(K_camera, R_ori_first, T_ori_first, P1_ori_first, P2_ori_first, Code_ori_first);
	for (int ii = 0; ii < xyz_ori_first.size(); ii++)
	{
		xyz_sfm_ori.push_back(xyz_ori_first[ii]);
	}
	if (xyz_sfm_ori.size() < 4)
	{
		oQProgressDialog.cancel();
		QMessageBox::warning(this, tr("SFM"), tr("Failed to initialize the posture in the predetermined direction: First KeyFrame points less than 4!"), QMessageBox::Ok, QMessageBox::Ok);
		cam_log_ui->hide();
		return;
	}

	Keypoints_in_queue.clear();
	Keypoints_contours_in_queue.clear();
	Keypoints_index_in_queue.clear();
	Keypoints_R_in_queue.clear();
	Keypoints_T_in_queue.clear();

	Keypoints_in_queue.push_back(KeyPoint_code_serial_copy[match_queue.top().position_x]);
	Keypoints_in_queue.push_back(KeyPoint_code_serial_copy[match_queue.top().position_y]);
	Keypoints_contours_in_queue.push_back(KeyPoint_contour_serial_copy[match_queue.top().position_x]);
	Keypoints_contours_in_queue.push_back(KeyPoint_contour_serial_copy[match_queue.top().position_y]);
	Keypoints_index_in_queue.push_back(KeyPoint_code_serial_index_copy[match_queue.top().position_x]);
	Keypoints_index_in_queue.push_back(KeyPoint_code_serial_index_copy[match_queue.top().position_y]);
	Keypoints_R_in_queue.push_back(cv::Mat::eye(3, 3, CV_64FC1));
	Keypoints_R_in_queue.push_back(R_ori_first);
	Keypoints_T_in_queue.push_back(cv::Mat::zeros(3, 1, CV_64FC1));
	Keypoints_T_in_queue.push_back(T_ori_first);

	remove_calculated_index_from_pair(KeyPoint_code_serial_copy, KeyPoint_contour_serial_copy, KeyPoint_code_serial_index_copy, match_queue.top());

	for (int pp = 0; pp < xyz_sfm_ori.size(); pp++)
	{
		xyz_sfm_ori[pp].weight = 0;
		for (int ss = 0; ss < Keypoints_in_queue.size(); ss++)
		{
			for (int tt = 0; tt < Keypoints_in_queue[ss].size(); tt++)
			{
				if (xyz_sfm_ori[pp].code == Keypoints_in_queue[ss][tt].code_num)
				{
					xyz_sfm_ori[pp].weight++;
					break;
				}
			}
		}
	}


	double* K = new double[5];
	double* Dis_K = new double[6];
	double* Dis_P = new double[2];
	double* Dis_T = new double[4];
	K[0] = fx;
	K[1] = fy;
	K[2] = cx;
	K[3] = cy;
	K[4] = 0;
	Dis_K[0] = 0;
	Dis_K[1] = 0;
	Dis_K[2] = 0;
	Dis_K[3] = 0;
	Dis_K[4] = 0;
	Dis_K[5] = 0;
	Dis_P[0] = 0;
	Dis_P[1] = 0;
	Dis_T[0] = 0;
	Dis_T[1] = 0;
	Dis_T[2] = 0;
	Dis_T[3] = 0;
	cv::Mat dis_coff = cv::Mat::eye(1, 12, CV_64FC1);
	bool need_optimize_C;
	bool need_optimize_F;
	double* range_f = new double[3];
	double* range_c = new double[3];
	if (SFM_setting_window->ui.unlimit_checkBox->isChecked())
	{
		delete[] range_f;
		delete[] range_c;
		range_f = nullptr;
		range_c = nullptr;
	}
	else
	{
		range_f[0] = fx;
		range_f[1] = fy;
		range_f[2] = ((double)SFM_setting_window->ui.limit_spinBox->value()) / 100.0 * fx < ((double)SFM_setting_window->ui.limit_spinBox->value()) / 100.0 * fy
			? ((double)SFM_setting_window->ui.limit_spinBox->value()) / 100.0 * fx : ((double)SFM_setting_window->ui.limit_spinBox->value()) / 100.0 * fy;
		range_c[0] = cx;
		range_c[1] = cy;
		range_c[2] = ((double)SFM_setting_window->ui.limit_spinBox->value()) / 100.0 * cx < ((double)SFM_setting_window->ui.limit_spinBox->value()) / 100.0 * cy
			? ((double)SFM_setting_window->ui.limit_spinBox->value()) / 100.0 * cx : ((double)SFM_setting_window->ui.limit_spinBox->value()) / 100.0 * cy;
	}


	while (KeyPoint_code_serial_copy.size() != 0)
	{
		oQProgressDialog.setValue(Keypoints_index_in_queue.size());
		std::priority_queue<match_3d_point_group> matched_3d;
		for (int ii = 0; ii < KeyPoint_code_serial_copy.size(); ii++)
		{
			int same_number = find_same_coded_point_to_3d(xyz_sfm_ori, KeyPoint_code_serial_copy[ii]);
			matched_3d.push(match_3d_point_group(same_number, ii));
		}
		std::vector<cv::Point3d> P_3d_temp;
		std::vector<cv::Point2d> P_ori_temp;
		std::vector<int> Code_ori_temp;
		cv::Mat R_ori_temp, T_ori_temp;
		cv::Mat mask_ori_temp;
		K_camera = cv::Mat(cv::Matx33d(
			K[0], K[4], K[2],
			0, K[1], K[3],
			0, 0, 1));
		dis_coff.at<double>(0) = Dis_K[0];
		dis_coff.at<double>(1) = Dis_K[1];
		dis_coff.at<double>(2) = Dis_P[0];
		dis_coff.at<double>(3) = Dis_P[1];
		dis_coff.at<double>(4) = Dis_K[2];
		dis_coff.at<double>(5) = Dis_K[3];
		dis_coff.at<double>(6) = Dis_K[4];
		dis_coff.at<double>(7) = Dis_K[5];
		dis_coff.at<double>(8) = Dis_T[0];
		dis_coff.at<double>(9) = Dis_T[2];
		dis_coff.at<double>(10) = Dis_T[1];
		dis_coff.at<double>(11) = Dis_T[3];
		bool debug = false;
		while (matched_3d.size())
		{
			if (matched_3d.top().matching_number < 4)
			{
				break;
			}
			if (cal_3d_transform(K_camera, dis_coff, xyz_sfm_ori, KeyPoint_code_serial_copy[matched_3d.top().position]
				, P_3d_temp, P_ori_temp, Code_ori_temp,R_ori_temp, T_ori_temp, mask_ori_temp) >= 3)
			{
				break;
			}
			matched_3d.pop();
		}
		if (matched_3d.size() == 0 || matched_3d.top().matching_number < 4)
		{
			break;
		}
		int new_n = 0;
		double min_ratio;
		update_points_in_point(Keypoints_in_queue,
			KeyPoint_code_serial_copy[matched_3d.top().position], K_camera, dis_coff,
			Keypoints_R_in_queue, Keypoints_T_in_queue, R_ori_temp, T_ori_temp,
			xyz_sfm_ori, min_ratio, new_n);
		Keypoints_in_queue.push_back(KeyPoint_code_serial_copy[matched_3d.top().position]);
		Keypoints_contours_in_queue.push_back(KeyPoint_contour_serial_copy[matched_3d.top().position]);
		Keypoints_index_in_queue.push_back(KeyPoint_code_serial_index_copy[matched_3d.top().position]);
		Keypoints_R_in_queue.push_back(R_ori_temp);
		Keypoints_T_in_queue.push_back(T_ori_temp);
		str_inf_for_log = 
			tr("Successfully add frame Image-") + QString::number(KeyPoint_code_serial_index_copy[matched_3d.top().position] + 1) + tr(": Minimum success rate of cross matching: ") +
			QString::number(min_ratio * 100) + tr("% ; Points added ") +
			QString::number(new_n) + " times." + tr(".\n") + str_inf_for_log;
		cam_log_ui->ui.textEdit->setText(str_inf_for_log);
		remove_calculated_index_from_index(KeyPoint_code_serial_copy, KeyPoint_contour_serial_copy, KeyPoint_code_serial_index_copy, matched_3d.top().position);
		//auto a = ((Keypoints_in_queue.size() - 2) % SFM_setting_window->ui.BA_interval_spinBox->value());
		for (int pp = 0; pp < xyz_sfm_ori.size(); pp++)
		{
			xyz_sfm_ori[pp].weight = 0;
		}
		for (int pp = 0; pp < xyz_sfm_ori.size(); pp++)
		{
			for (int ss = 0; ss < Keypoints_in_queue.size(); ss++)
			{
				for (int tt = 0; tt < Keypoints_in_queue[ss].size(); tt++)
				{
					if (xyz_sfm_ori[pp].code == Keypoints_in_queue[ss][tt].code_num)
					{
						xyz_sfm_ori[pp].weight++;
						break;
					}
				}
			}
		}

		//std::vector<Eigen::Quaterniond> camQvec_part;
		//std::vector<Eigen::Vector3d> camTvec_part;
		//std::vector<std::vector<cv::Point2d>> Re_project_Map_part;
		//for (int pp = 0; pp < Keypoints_R_in_queue.size(); pp++)
		//{
		//	Eigen::Matrix<double, 3, 3> matrix_temp;
		//	cv::cv2eigen(Keypoints_R_in_queue[pp], matrix_temp);
		//	camQvec_part.push_back(Eigen::Quaterniond(matrix_temp));
		//	Eigen::Matrix<double, 3, 1> matrix_temp2;
		//	cv::cv2eigen(Keypoints_T_in_queue[pp], matrix_temp2);
		//	camTvec_part.push_back(Eigen::Vector3d(matrix_temp2));
		//}
		//need_optimize_C = true;
		//if (SFM_setting_window->ui.fixed_c_checkBox->isChecked())
		//{
		//	need_optimize_C = false;
		//}
		//else if (SFM_setting_window->ui.optimize_checkBox->isChecked())
		//{
		//	need_optimize_C = false;
		//}
		//else
		//{
		//	if (Keypoints_in_queue.size() < SFM_setting_window->ui.min_cxcy_BA_spinBox->value())
		//	{
		//		need_optimize_C = false;
		//	}
		//}
		//need_optimize_F = true;
		//if (SFM_setting_window->ui.fixed_f_checkBox->isChecked())
		//{
		//	need_optimize_F = false;
		//}
		//else if (SFM_setting_window->ui.optimize_checkBox->isChecked())
		//{
		//	need_optimize_F = false;
		//}
		//else
		//{
		//	if (Keypoints_in_queue.size() < SFM_setting_window->ui.min_cxcy_BA_spinBox->value())
		//	{
		//		need_optimize_F = false;
		//	}
		//}
		//ceres::Solver::Summary summary_part;
		//optimize_SFM_incremental_thread thread_cal_part(Keypoints_in_queue, xyz_sfm_ori, camQvec_part, camTvec_part, Re_project_Map_part, K, Dis_K, Dis_P, Dis_T, &summary_part
		//	, SFM_setting_window->ui.max_iter_spinBox->value()
		//	, ((double)SFM_setting_window->ui.stop_spinBox->value())* pow(10, SFM_setting_window->ui.stop_spinBox_2->value())
		//	, SFM_setting_window->ui.thread_spinBox->value()
		//	, SFM_setting_window->ui.timeout_spinBox->value()
		//	, SFM_setting_window->ui.loss_comboBox->currentIndex(), SFM_setting_window->ui.loss_doubleSpinBox->value()
		//	, SFM_setting_window->ui.k_distort_comboBox->currentText().toInt()
		//	, SFM_setting_window->ui.p_distort_comboBox->currentText().toInt()
		//	, SFM_setting_window->ui.s_distort_comboBox->currentText().toInt()
		//	, SFM_setting_window->ui.fixed_camera_checkBox->isChecked()
		//	, need_optimize_F
		//	, need_optimize_C
		//	, SFM_setting_window->ui.skew_checkBox->isChecked()
		//	, SFM_setting_window->ui.uniform_f_checkBox->isChecked(), range_f, range_c, true);
		//thread_cal_part.start();
		//while (1)
		//{
		//	if (!thread_cal_part.isRunning())
		//	{
		//		break;
		//	}
		//	QString inf_log = "";
		//	for (int ss = summary_part.iterations.size() - 1; ss >= 0; ss--)
		//	{
		//		inf_log = inf_log + QString("%1").arg(ss + 1, 4, 10, QLatin1Char('0')) + tr(". ");
		//		inf_log = inf_log + tr("[Cost]:") + QString::number(summary_part.iterations[ss].cost, 'e', 3) + tr("; ");
		//		inf_log = inf_log + tr("[Gradient Norm]:") + QString::number(summary_part.iterations[ss].gradient_norm, 'e', 3) + tr("; ");
		//		inf_log = inf_log + tr("[Step Norm]:") + QString::number(summary_part.iterations[ss].step_norm, 'e', 3) + tr("; ");
		//		inf_log = inf_log + tr("[Iteration Time(second)]:") + QString::number(summary_part.iterations[ss].iteration_time_in_seconds, 'e', 3) + tr(";");
		//		inf_log = inf_log + tr("\n");
		//	}
		//	cam_log_ui->ui.textEdit->setText(inf_log);
		//	loop_sleep(20);
		//}
		//if (KeyPoint_code_serial_index_copy.size()==0)
		//{
		//	int aaaa = 1;
		//}
		//xyz_sfm_ori = thread_cal_part.point_clouds_thread;
		//camQvec_part = thread_cal_part.camQvec_thread;
		//camTvec_part = thread_cal_part.camTvec_thread;
		//Re_project_Map_part = thread_cal_part.Re_project_Map_thread;
		//K = thread_cal_part.camK_thread;
		//Dis_K = thread_cal_part.camDK_thread;
		//Dis_P = thread_cal_part.camDP_thread;
		//Dis_T = thread_cal_part.camDT_thread;

		//for (int pp = 1; pp < Keypoints_R_in_queue.size(); pp++)
		//{
		//	cv::Mat mat_temp;
		//	cv::eigen2cv(camQvec_part[pp].matrix(), mat_temp);
		//	Keypoints_R_in_queue[pp] = mat_temp;
		//	cv::Mat mat_temp2;
		//	cv::eigen2cv(camTvec_part[pp], mat_temp2);
		//	Keypoints_T_in_queue[pp] = mat_temp2;
		//}

		//str_inf_for_log =
		//	tr("Optimized partly camera parameters:") +
		//	tr("Fx=") + QString::number(K[0]) + tr("; ") +
		//	tr("Fy=") + QString::number(K[1]) + tr("; ") +
		//	tr("Cx=") + QString::number(K[2]) + tr("; ") +
		//	tr("Cy=") + QString::number(K[3]) + tr("; ") +
		//	tr("Skew=") + QString::number(K[4]) + tr("; ") +
		//	tr("Distortion_K=[") + QString::number(Dis_K[0]) + tr(",") + QString::number(Dis_K[1]) + tr(",") + QString::number(Dis_K[2]) + tr(",") + QString::number(Dis_K[3]) + tr(",") +
		//	QString::number(Dis_K[4]) + tr(",") + QString::number(Dis_K[5]) + tr("]; ") +
		//	tr("Distortion_P=[") + QString::number(Dis_P[0]) + tr(",") + QString::number(Dis_P[1]) + tr("]; ") +
		//	tr("Distortion_T=[") + QString::number(Dis_T[0]) + tr(",") + QString::number(Dis_T[1]) + tr(",") + QString::number(Dis_T[2]) + tr(",") + QString::number(Dis_T[3]) + tr("].\n")
		//	+ str_inf_for_log;
		//cam_log_ui->ui.textEdit->setText(str_inf_for_log);
		if (((Keypoints_in_queue.size() - 2) % SFM_setting_window->ui.BA_interval_spinBox->value()) == 0 && xyz_sfm_ori.size() != 0)
		{
			std::vector<Eigen::Quaterniond> camQvec;
			std::vector<Eigen::Vector3d> camTvec;
			std::vector<std::vector<cv::Point2d>> Re_project_Map;
			for (int pp = 0; pp < Keypoints_R_in_queue.size(); pp++)
			{
				Eigen::Matrix<double, 3, 3> matrix_temp;
				cv::cv2eigen(Keypoints_R_in_queue[pp], matrix_temp);
				camQvec.push_back(Eigen::Quaterniond(matrix_temp));
				camQvec[camQvec.size() - 1].normalize();
				Eigen::Matrix<double, 3, 1> matrix_temp2;
				cv::cv2eigen(Keypoints_T_in_queue[pp], matrix_temp2);
				camTvec.push_back(Eigen::Vector3d(matrix_temp2));
			}
			need_optimize_C = true;
			if (SFM_setting_window->ui.fixed_c_checkBox->isChecked())
			{
				need_optimize_C = false;
			}
			else if (SFM_setting_window->ui.optimize_checkBox->isChecked())
			{
				need_optimize_C = false;
			}
			else
			{
				if (Keypoints_in_queue.size() < SFM_setting_window->ui.min_cxcy_BA_spinBox->value())
				{
					need_optimize_C = false;
				}
			}
			need_optimize_F = true;
			if (SFM_setting_window->ui.fixed_f_checkBox->isChecked())
			{
				need_optimize_F = false;
			}
			else if (SFM_setting_window->ui.optimize_checkBox->isChecked())
			{
				need_optimize_F = false;
			}
			else
			{
				if (Keypoints_in_queue.size() < SFM_setting_window->ui.min_cxcy_BA_spinBox->value())
				{
					need_optimize_F = false;
				}
			}
			ceres::Solver::Summary summary;
			optimize_SFM_incremental_thread thread_cal(Keypoints_in_queue, xyz_sfm_ori, camQvec, camTvec, Re_project_Map, K, Dis_K, Dis_P, Dis_T, &summary
				, SFM_setting_window->ui.max_iter_spinBox->value()
				, ((double)SFM_setting_window->ui.stop_spinBox->value())* pow(10, SFM_setting_window->ui.stop_spinBox_2->value())
				, SFM_setting_window->ui.thread_spinBox->value()
				, SFM_setting_window->ui.timeout_spinBox->value()
				, SFM_setting_window->ui.loss_comboBox->currentIndex(), SFM_setting_window->ui.loss_doubleSpinBox->value()
				, SFM_setting_window->ui.k_distort_comboBox->currentText().toInt()
				, SFM_setting_window->ui.p_distort_comboBox->currentText().toInt()
				, SFM_setting_window->ui.s_distort_comboBox->currentText().toInt()
				, SFM_setting_window->ui.fixed_camera_checkBox->isChecked()
				, need_optimize_F
				, need_optimize_C
				, SFM_setting_window->ui.skew_checkBox->isChecked()
				, SFM_setting_window->ui.uniform_f_checkBox->isChecked(), range_f, range_c);
			thread_cal.start();
			while (1)
			{
				if (!thread_cal.isRunning())
				{
					break;
				}
				QString inf_log = "";
				for (int ss = summary.iterations.size() - 1; ss >= 0; ss--)
				{
					inf_log = inf_log + QString("%1").arg(ss + 1, 4, 10, QLatin1Char('0')) + tr(". ");
					inf_log = inf_log + tr("[Cost]:") + QString::number(summary.iterations[ss].cost, 'e', 3) + tr("; ");
					inf_log = inf_log + tr("[Gradient Norm]:") + QString::number(summary.iterations[ss].gradient_norm, 'e', 3) + tr("; ");
					inf_log = inf_log + tr("[Step Norm]:") + QString::number(summary.iterations[ss].step_norm, 'e', 3) + tr("; ");
					inf_log = inf_log + tr("[Iteration Time(second)]:") + QString::number(summary.iterations[ss].iteration_time_in_seconds, 'e', 3) + tr(";");
					inf_log = inf_log + tr("\n");
				}
				cam_log_ui->ui.textEdit->setText(inf_log);
				loop_sleep(20);
			}
			xyz_sfm_ori = thread_cal.point_clouds_thread;
			camQvec = thread_cal.camQvec_thread;
			camTvec = thread_cal.camTvec_thread;
			Re_project_Map = thread_cal.Re_project_Map_thread;
			K = thread_cal.camK_thread;
			Dis_K = thread_cal.camDK_thread;
			Dis_P = thread_cal.camDP_thread;
			Dis_T = thread_cal.camDT_thread;

			for (int pp = 1; pp < Keypoints_R_in_queue.size(); pp++)
			{
				cv::Mat mat_temp;
				cv::eigen2cv(camQvec[pp].matrix(), mat_temp);
				Keypoints_R_in_queue[pp] = mat_temp;
				cv::Mat mat_temp2;
				cv::eigen2cv(camTvec[pp], mat_temp2);
				Keypoints_T_in_queue[pp] = mat_temp2;
			}

			str_inf_for_log =
				tr("Optimized camera parameters:") + 
				tr("Fx=") + QString::number(K[0]) + tr("; ") +
				tr("Fy=") + QString::number(K[1]) + tr("; ") +
				tr("Cx=") + QString::number(K[2]) + tr("; ") +
				tr("Cy=") + QString::number(K[3]) + tr("; ") +
				tr("Skew=") + QString::number(K[4]) + tr("; ") +
				tr("Distortion_K=[") + QString::number(Dis_K[0]) + tr(",") + QString::number(Dis_K[1]) + tr(",") + QString::number(Dis_K[2]) + tr(",") + QString::number(Dis_K[3]) + tr(",") +
				QString::number(Dis_K[4]) + tr(",") + QString::number(Dis_K[5]) + tr("]; ") +
				tr("Distortion_P=[") + QString::number(Dis_P[0]) + tr(",") + QString::number(Dis_P[1]) + tr("]; ") +
				tr("Distortion_T=[") + QString::number(Dis_T[0]) + tr(",") + QString::number(Dis_T[1]) + tr(",") + QString::number(Dis_T[2]) + tr(",") + QString::number(Dis_T[3]) + tr("].\n")
				+ str_inf_for_log;
			cam_log_ui->ui.textEdit->setText(str_inf_for_log);
		}

		Result_xyz_serial = xyz_sfm_ori;
		xyz_need_levelling_serial.clear();
		xyz_need_levelling_code_serial.clear();
		for (int qs = 0; qs < Result_xyz_serial.size(); qs++)
		{
			xyz_need_levelling_serial.push_back(false);
			xyz_need_levelling_code_serial.push_back(Result_xyz_serial[qs].code);
		}
		std::sort(xyz_need_levelling_code_serial.begin(), xyz_need_levelling_code_serial.end());
		update_points_view();
		qApp->processEvents();
	}
	oQProgressDialog.setValue(Keypoints_index_in_queue.size() + 1);

	//sbsbsb1
	/*std::vector<cv::Point3d> out_3d;
	std::vector<cv::Point2d> out_2d;
	cv::Mat R_temp, T_temp, mask_temp;
	for (int ii = 0; ii < 5; ii++)
	{
		for (int jj = 0; jj < 5; jj++)
		{
			xyz_sfm_ori.push_back(sfm_3d_group(ii + jj * 5 + 1, 100, ii * 250, jj * 250, 0, 0));
		}
	}
	out_3d.clear();
	out_2d.clear();
	for (int ii = 0; ii < xyz_sfm_ori.size(); ii++)
	{
		for (int jj = 0; jj < KeyPoint_code_serial_copy[0].size(); jj++)
		{
			if (xyz_sfm_ori[ii].code == KeyPoint_code_serial_copy[0][jj].code_num)
			{
				out_3d.push_back(cv::Point3d(xyz_sfm_ori[ii].x, xyz_sfm_ori[ii].y, xyz_sfm_ori[ii].z));
				out_2d.push_back(cv::Point2d(KeyPoint_code_serial_copy[0][jj].x, KeyPoint_code_serial_copy[0][jj].y));
			}
		}
	}
	cv::solvePnPRansac(out_3d, out_2d, K_camera, dis_coff, R_temp, T_temp, false, 100, 8.0, 0.99, mask_temp);
	Rodrigues(R_temp, R_temp);
	if ((double)countNonZero(mask_temp)  < 20)
	{
		return;
	}
	for (int ii = 0; ii < xyz_sfm_ori.size(); ii++)
	{
		double x = xyz_sfm_ori[ii].x;
		double y = xyz_sfm_ori[ii].y;
		double z = xyz_sfm_ori[ii].z;
		xyz_sfm_ori[ii].x = R_temp.at<double>(0, 0) * x + R_temp.at<double>(0, 1) * y + R_temp.at<double>(0, 2) * z + T_temp.at<double>(0, 0);
		xyz_sfm_ori[ii].y = R_temp.at<double>(1, 0) * x + R_temp.at<double>(1, 1) * y + R_temp.at<double>(1, 2) * z + T_temp.at<double>(1, 0);
		xyz_sfm_ori[ii].z = R_temp.at<double>(2, 0) * x + R_temp.at<double>(2, 1) * y + R_temp.at<double>(2, 2) * z + T_temp.at<double>(2, 0);
	}

	out_3d.clear();
	out_2d.clear();
	for (int ii = 0; ii < xyz_sfm_ori.size(); ii++)
	{
		for (int jj = 0; jj < KeyPoint_code_serial_copy[0].size(); jj++)
		{
			if (xyz_sfm_ori[ii].code == KeyPoint_code_serial_copy[0][jj].code_num)
			{
				out_3d.push_back(cv::Point3d(xyz_sfm_ori[ii].x, xyz_sfm_ori[ii].y, xyz_sfm_ori[ii].z));
				out_2d.push_back(cv::Point2d(KeyPoint_code_serial_copy[0][jj].x, KeyPoint_code_serial_copy[0][jj].y));
			}
		}
	}
	cv::solvePnPRansac(out_3d, out_2d, K_camera, dis_coff, R_temp, T_temp, false, 100, 8.0, 0.99, mask_temp);
	Rodrigues(R_temp, R_temp);

	Keypoints_in_queue.push_back(KeyPoint_code_serial_copy[0]);
	Keypoints_contours_in_queue.push_back(KeyPoint_contour_serial_copy[0]);
	Keypoints_index_in_queue.push_back(KeyPoint_code_serial_index_copy[0]);
	Keypoints_R_in_queue.push_back(R_temp);
	Keypoints_T_in_queue.push_back(T_temp);


	out_3d.clear();
	out_2d.clear();
	for (int ii = 0; ii < xyz_sfm_ori.size(); ii++)
	{
		for (int jj = 0; jj < KeyPoint_code_serial_copy[1].size(); jj++)
		{
			if (xyz_sfm_ori[ii].code == KeyPoint_code_serial_copy[1][jj].code_num)
			{
				out_3d.push_back(cv::Point3d(xyz_sfm_ori[ii].x, xyz_sfm_ori[ii].y, xyz_sfm_ori[ii].z));
				out_2d.push_back(cv::Point2d(KeyPoint_code_serial_copy[1][jj].x, KeyPoint_code_serial_copy[1][jj].y));
			}
		}
	}
	cv::solvePnPRansac(out_3d, out_2d, K_camera, dis_coff, R_temp, T_temp, false, 100, 8.0, 0.99, mask_temp);
	Rodrigues(R_temp, R_temp);
	double module_val = sqrt(T_temp.at<double>(0, 0) * T_temp.at<double>(0, 0) + T_temp.at<double>(1, 0) * T_temp.at<double>(1, 0) + T_temp.at<double>(2, 0) * T_temp.at<double>(2, 0));
	for (int ii = 0; ii < xyz_sfm_ori.size(); ii++)
	{
		xyz_sfm_ori[ii].x /= module_val;
		xyz_sfm_ori[ii].y /= module_val;
		xyz_sfm_ori[ii].z /= module_val;
	}
	out_3d.clear();
	out_2d.clear();
	for (int ii = 0; ii < xyz_sfm_ori.size(); ii++)
	{
		for (int jj = 0; jj < KeyPoint_code_serial_copy[1].size(); jj++)
		{
			if (xyz_sfm_ori[ii].code == KeyPoint_code_serial_copy[1][jj].code_num)
			{
				out_3d.push_back(cv::Point3d(xyz_sfm_ori[ii].x, xyz_sfm_ori[ii].y, xyz_sfm_ori[ii].z));
				out_2d.push_back(cv::Point2d(KeyPoint_code_serial_copy[1][jj].x, KeyPoint_code_serial_copy[1][jj].y));
			}
		}
	}
	cv::solvePnPRansac(out_3d, out_2d, K_camera, dis_coff, R_temp, T_temp, false, 100, 8.0, 0.99, mask_temp);
	Rodrigues(R_temp, R_temp);
	Keypoints_in_queue.push_back(KeyPoint_code_serial_copy[1]);
	Keypoints_contours_in_queue.push_back(KeyPoint_contour_serial_copy[1]);
	Keypoints_index_in_queue.push_back(KeyPoint_code_serial_index_copy[1]);
	Keypoints_R_in_queue.push_back(R_temp);
	Keypoints_T_in_queue.push_back(T_temp);

	for (int pp = 2; pp < KeyPoint_code_serial_copy.size(); pp++)
	{
		out_3d.clear();
		out_2d.clear();
		for (int ii = 0; ii < xyz_sfm_ori.size(); ii++)
		{
			for (int jj = 0; jj < KeyPoint_code_serial_copy[pp].size(); jj++)
			{
				if (xyz_sfm_ori[ii].code == KeyPoint_code_serial_copy[pp][jj].code_num)
				{
					out_3d.push_back(cv::Point3d(xyz_sfm_ori[ii].x, xyz_sfm_ori[ii].y, xyz_sfm_ori[ii].z));
					out_2d.push_back(cv::Point2d(KeyPoint_code_serial_copy[pp][jj].x, KeyPoint_code_serial_copy[pp][jj].y));
				}
			}
		}
		cv::solvePnPRansac(out_3d, out_2d, K_camera, dis_coff, R_temp, T_temp, false, 100, 8.0, 0.99, mask_temp);
		Rodrigues(R_temp, R_temp);
		if ((double)countNonZero(mask_temp) < 20)
		{
			continue;
		}
		Keypoints_in_queue.push_back(KeyPoint_code_serial_copy[pp]);
		Keypoints_contours_in_queue.push_back(KeyPoint_contour_serial_copy[pp]);
		Keypoints_index_in_queue.push_back(KeyPoint_code_serial_index_copy[pp]);
		Keypoints_R_in_queue.push_back(R_temp);
		Keypoints_T_in_queue.push_back(T_temp);
	}
	std::vector<Eigen::Quaterniond> camQvec;
	std::vector<Eigen::Vector3d> camTvec;
	std::vector<std::vector<cv::Point2d>> Re_project_Map;
	for (int pp = 0; pp < Keypoints_R_in_queue.size(); pp++)
	{
		Eigen::Matrix<double, 3, 3> matrix_temp;
		cv::cv2eigen(Keypoints_R_in_queue[pp], matrix_temp);
		camQvec.push_back(Eigen::Quaterniond(matrix_temp));
		Eigen::Matrix<double, 3, 1> matrix_temp2;
		cv::cv2eigen(Keypoints_T_in_queue[pp], matrix_temp2);
		camTvec.push_back(Eigen::Vector3d(matrix_temp2));
	}

	std::vector<sfm_3d_group> xyz_sfm_ori_copy = xyz_sfm_ori;
	xyz_sfm_ori.clear();
	for (int ii = 0; ii < xyz_sfm_ori_copy.size(); ii++)
	{
		int flag = 0;
		for (int pp = 0; pp < KeyPoint_code_serial_copy.size(); pp++)
		{
			for (int qq = 0; qq < KeyPoint_code_serial_copy[pp].size(); qq++)
			{
				if (KeyPoint_code_serial_copy[pp][qq].code_num == xyz_sfm_ori_copy[ii].code)
				{
					flag++;
				}
			}
		}
		if (flag>=2)
		{
			xyz_sfm_ori.push_back(xyz_sfm_ori_copy[ii]);
		}
	}*/

    //sbsbsb2
	/*Keypoints_R_in_queue.push_back(R_save[index_for_L]);
	Keypoints_R_in_queue.push_back(R_save[index_for_R]);
	Keypoints_T_in_queue.push_back(T_save[index_for_L]);
	Keypoints_T_in_queue.push_back(T_save[index_for_R]);
	Keypoints_in_queue.push_back(KeyPoint_code_serial_copy[0]);
	Keypoints_contours_in_queue.push_back(KeyPoint_contour_serial_copy[0]);
	Keypoints_index_in_queue.push_back(KeyPoint_code_serial_index_copy[0]);
	Keypoints_in_queue.push_back(KeyPoint_code_serial_copy[1]);
	Keypoints_contours_in_queue.push_back(KeyPoint_contour_serial_copy[1]);
	Keypoints_index_in_queue.push_back(KeyPoint_code_serial_index_copy[1]);
	std::vector<Eigen::Quaterniond> camQvec;
	std::vector<Eigen::Vector3d> camTvec;
	std::vector<std::vector<cv::Point2d>> Re_project_Map;
	for (int pp = 0; pp < Keypoints_R_in_queue.size(); pp++)
	{
		Eigen::Matrix<double, 3, 3> matrix_temp;
		cv::cv2eigen(Keypoints_R_in_queue[pp], matrix_temp);
		camQvec.push_back(Eigen::Quaterniond(matrix_temp));
		Eigen::Matrix<double, 3, 1> matrix_temp2;
		cv::cv2eigen(Keypoints_T_in_queue[pp], matrix_temp2);
		camTvec.push_back(Eigen::Vector3d(matrix_temp2));
	}*/



	//sbsbsb
	std::vector<Eigen::Quaterniond> camQvec;
	std::vector<Eigen::Vector3d> camTvec;
	std::vector<std::vector<cv::Point2d>> Re_project_Map;
	for (int pp = 0; pp < Keypoints_R_in_queue.size(); pp++)
	{
		Eigen::Matrix<double, 3, 3> matrix_temp;
		cv::cv2eigen(Keypoints_R_in_queue[pp], matrix_temp);
		camQvec.push_back(Eigen::Quaterniond(matrix_temp));
		camQvec[camQvec.size() - 1].normalize();
		Eigen::Matrix<double, 3, 1> matrix_temp2;
		cv::cv2eigen(Keypoints_T_in_queue[pp], matrix_temp2);
		camTvec.push_back(Eigen::Vector3d(matrix_temp2));
	}
	std::vector<sfm_3d_group> xyz_sfm_ori_copy = xyz_sfm_ori;
	xyz_sfm_ori.clear();
	for (int ii = 0; ii < xyz_sfm_ori_copy.size(); ii++)
	{
		if (xyz_sfm_ori_copy[ii].weight >= SFM_setting_window->ui.min_view_spinBox->value())
		{
			xyz_sfm_ori.push_back(xyz_sfm_ori_copy[ii]);
		}
	}

	if (xyz_sfm_ori.size() == 0)
	{
		oQProgressDialog.cancel();
		QMessageBox::warning(this, tr("SFM"), tr("No points can be reconstructed or pass [Minimum-View] threshold!"), QMessageBox::Ok, QMessageBox::Ok);
		
		Result_xyz_serial.clear();
		Result_xyz_after_op_serial.clear();
		Result_R_after_op.clear();
		Result_T_after_op.clear();
		Result_Ori_after_op.clear();
		xyz_need_levelling_serial.clear();
		xyz_need_levelling_code_serial.clear();
		Result_R.clear();
		Result_T.clear();
		Result_Ori.clear();
		Result_err_every.clear();
		Keypoints_in_queue.clear();
		Keypoints_contours_in_queue.clear();
		Keypoints_index_in_queue.clear();
		Keypoints_R_in_queue.clear();
		Keypoints_T_in_queue.clear();
		update_points_view();
		cam_log_ui->hide();
		return;
	}
	need_optimize_C = true;
	if (SFM_setting_window->ui.fixed_c_checkBox->isChecked())
	{
		need_optimize_C = false;
	}
	need_optimize_F = true;
	if (SFM_setting_window->ui.fixed_f_checkBox->isChecked())
	{
		need_optimize_F = false;
	}
	ceres::Solver::Summary summary;
	optimize_SFM_incremental_thread thread_cal(Keypoints_in_queue, xyz_sfm_ori, camQvec, camTvec, Re_project_Map, K, Dis_K, Dis_P, Dis_T, &summary
		, SFM_setting_window->ui.max_iter_spinBox->value()
		, ((double)SFM_setting_window->ui.stop_spinBox->value())* pow(10, SFM_setting_window->ui.stop_spinBox_2->value())
		, SFM_setting_window->ui.thread_spinBox->value()
		, SFM_setting_window->ui.timeout_spinBox->value()
		, SFM_setting_window->ui.loss_comboBox->currentIndex(), SFM_setting_window->ui.loss_doubleSpinBox->value()
		, SFM_setting_window->ui.k_distort_comboBox->currentText().toInt()
		, SFM_setting_window->ui.p_distort_comboBox->currentText().toInt()
		, SFM_setting_window->ui.s_distort_comboBox->currentText().toInt()
		, SFM_setting_window->ui.fixed_camera_checkBox->isChecked()
		, need_optimize_F
		, need_optimize_C
		, SFM_setting_window->ui.skew_checkBox->isChecked()
		, SFM_setting_window->ui.uniform_f_checkBox->isChecked(), range_f, range_c);
	thread_cal.start();
	while (1)
	{
		if (!thread_cal.isRunning())
		{
			break;
		}
		QString inf_log = "";
		for (int ss = summary.iterations.size() - 1; ss >= 0; ss--)
		{
			inf_log = inf_log + QString("%1").arg(ss + 1, 4, 10, QLatin1Char('0')) + tr(". ");
			inf_log = inf_log + tr("[Cost]:") + QString::number(summary.iterations[ss].cost, 'e', 3) + tr("; ");
			inf_log = inf_log + tr("[Gradient Norm]:") + QString::number(summary.iterations[ss].gradient_norm, 'e', 3) + tr("; ");
			inf_log = inf_log + tr("[Step Norm]:") + QString::number(summary.iterations[ss].step_norm, 'e', 3) + tr("; ");
			inf_log = inf_log + tr("[Iteration Time(second)]:") + QString::number(summary.iterations[ss].iteration_time_in_seconds, 'e', 3) + tr(";");
			inf_log = inf_log + tr("\n");
		}
		cam_log_ui->ui.textEdit->setText(inf_log);
		loop_sleep(20);
	}

	xyz_sfm_ori = thread_cal.point_clouds_thread;
	camQvec = thread_cal.camQvec_thread;
	camTvec = thread_cal.camTvec_thread;
	Re_project_Map = thread_cal.Re_project_Map_thread;
	K = thread_cal.camK_thread;
	Dis_K = thread_cal.camDK_thread;
	Dis_P = thread_cal.camDP_thread;
	Dis_T = thread_cal.camDT_thread;
	str_inf_for_log =
		tr("Global BA completion.\n") + str_inf_for_log;
	cam_log_ui->ui.textEdit->setText(str_inf_for_log);
	Result_Camera_K[0] = K[0];
	Result_Camera_K[1] = K[1];
	Result_Camera_K[2] = K[2];
	Result_Camera_K[3] = K[3];
	Result_Camera_K[4] = K[4];
	Result_Camera_Dis_K[0] = Dis_K[0];
	Result_Camera_Dis_K[1] = Dis_K[1];
	Result_Camera_Dis_K[2] = Dis_K[2];
	Result_Camera_Dis_K[3] = Dis_K[3];
	Result_Camera_Dis_K[4] = Dis_K[4];
	Result_Camera_Dis_K[5] = Dis_K[5];
	Result_Camera_Dis_P[0] = Dis_P[0];
	Result_Camera_Dis_P[1] = Dis_P[1];
	Result_Camera_Dis_T[0] = Dis_T[0];
	Result_Camera_Dis_T[1] = Dis_T[1];
	Result_Camera_Dis_T[2] = Dis_T[2];
	Result_Camera_Dis_T[3] = Dis_T[3];
	for (int pp = 1; pp < Keypoints_R_in_queue.size(); pp++)
	{
		cv::Mat mat_temp;
		cv::eigen2cv(camQvec[pp].matrix(), mat_temp);
		Keypoints_R_in_queue[pp] = mat_temp;
		cv::Mat mat_temp2;
		cv::eigen2cv(camTvec[pp], mat_temp2);
		Keypoints_T_in_queue[pp] = mat_temp2;
	}
	Result_err_every.clear();
	Result_R.clear();
	Result_T.clear();
	for (int ii = 0; ii < KeyPoint_code_serial.size(); ii++)
	{
		bool exist = false;
		for (int jj = 0; jj < Keypoints_index_in_queue.size(); jj++)
		{
			if (ii == Keypoints_index_in_queue[jj])
			{
				exist = true;
				Result_err_every.push_back(Re_project_Map[jj]);
				Result_R.push_back(camQvec[jj]);
				Result_T.push_back(camTvec[jj]);
				break;
			}
		}
		if (!exist)
		{
			std::vector<cv::Point2d> temp_err = {};
			for (int jj = 0; jj < KeyPoint_code_serial[ii].size(); jj++)
			{
				temp_err.push_back(cv::Point2d(std::numeric_limits<double>::max(), std::numeric_limits<double>::max()));
			}
			Result_err_every.push_back(temp_err);
			Result_R.push_back(Eigen::Quaterniond(0, 0, 0, 0));
			Result_T.push_back(Eigen::Vector3d(0, 0, 0));
		}
	}


	K_camera = cv::Mat(cv::Matx33d(
		K[0], K[4], K[2],
		0, K[1], K[3],
		0, 0, 1));
	dis_coff.at<double>(0) = Dis_K[0];
	dis_coff.at<double>(1) = Dis_K[1];
	dis_coff.at<double>(2) = Dis_P[0];
	dis_coff.at<double>(3) = Dis_P[1];
	dis_coff.at<double>(4) = Dis_K[2];
	dis_coff.at<double>(5) = Dis_K[3];
	dis_coff.at<double>(6) = Dis_K[4];
	dis_coff.at<double>(7) = Dis_K[5];
	dis_coff.at<double>(8) = Dis_T[0];
	dis_coff.at<double>(9) = Dis_T[2];
	dis_coff.at<double>(10) = Dis_T[1];
	dis_coff.at<double>(11) = Dis_T[3];

	for (int ii = 0; ii < Result_contours_serial.size();ii++)
	{
		for (int jj = 0; jj < Result_contours_serial[ii].size(); jj++)
		{
			ellipse_view_pose pose_now = calculate_ellipse(Result_contours_serial[ii][jj], K_camera, dis_coff, Result_R[ii], Result_T[ii]);
			KeyPoint_code_serial[ii][jj].pose = pose_now;
		}
	}
	for (int ii = 0; ii < xyz_sfm_ori.size(); ii++)
	{
		std::vector<Eigen::Vector3d> vec_3d_temp;
		for (int pp = 0; pp < KeyPoint_code_serial.size(); pp++)
		{
			for (int qq = 0; qq < KeyPoint_code_serial[pp].size(); qq++)
			{
				if (KeyPoint_code_serial[pp][qq].code_num == xyz_sfm_ori[ii].code)
				{

					Eigen::Vector3d camTvec_1;
					camTvec_1.x() = KeyPoint_code_serial[pp][qq].pose.view_N_1[0];
					camTvec_1.y() = KeyPoint_code_serial[pp][qq].pose.view_N_1[1];
					camTvec_1.z() = KeyPoint_code_serial[pp][qq].pose.view_N_1[2];
					Eigen::Vector3d camTvec_2;
					camTvec_2.x() = KeyPoint_code_serial[pp][qq].pose.view_N_2[0];
					camTvec_2.y() = KeyPoint_code_serial[pp][qq].pose.view_N_2[1];
					camTvec_2.z() = KeyPoint_code_serial[pp][qq].pose.view_N_2[2];


					Eigen::Vector3d Cres_1 = camTvec_1.matrix().transpose() * Result_R[pp].matrix();
					Eigen::Vector3d Cres_2 = camTvec_2.matrix().transpose() * Result_R[pp].matrix();
					vec_3d_temp.push_back(Cres_1);
					vec_3d_temp.push_back(Cres_2);
					break;
				}
			}
		}
		double flag = std::numeric_limits<double>::max();
		Eigen::Vector3d vec_now;
		for (int pp = 0; pp < vec_3d_temp.size(); pp++)
		{
			double val_now = 0;
			for (int qq = 0; qq < (vec_3d_temp.size() / 2); qq++)
			{
				if (pp / 2 == qq)
				{
					continue;
				}
				double f_1 = abs(vec_3d_temp[qq * 2].x() - vec_3d_temp[pp].x()) + abs(vec_3d_temp[qq * 2].y() - vec_3d_temp[pp].y()) + abs(vec_3d_temp[qq * 2].z() - vec_3d_temp[pp].z());
				double f_2 = abs(vec_3d_temp[qq * 2 + 1].x() - vec_3d_temp[pp].x()) + abs(vec_3d_temp[qq * 2 + 1].y() - vec_3d_temp[pp].y()) + abs(vec_3d_temp[qq * 2 + 1].z() - vec_3d_temp[pp].z());
				double f_3 = abs(vec_3d_temp[qq * 2].x() + vec_3d_temp[pp].x()) + abs(vec_3d_temp[qq * 2].y() + vec_3d_temp[pp].y()) + abs(vec_3d_temp[qq * 2].z() + vec_3d_temp[pp].z());
				double f_4 = abs(vec_3d_temp[qq * 2 + 1].x() + vec_3d_temp[pp].x()) + abs(vec_3d_temp[qq * 2 + 1].y() + vec_3d_temp[pp].y()) + abs(vec_3d_temp[qq * 2 + 1].z() + vec_3d_temp[pp].z());
				double f = f_1;
				f = f > f_2 ? f_2 : f;
				f = f > f_3 ? f_3 : f;
				f = f > f_4 ? f_4 : f;
				val_now += f;
			}
			if (val_now < flag)
			{
				flag = val_now;
				vec_now = vec_3d_temp[pp];
			}
		}

		Result_Ori.push_back(vec_now);
	}
	KeyPoint_code_serial_update = KeyPoint_code_serial;
	if (SFM_setting_window->ui.circle_correctin_checkBox->isChecked())
	{
		for (int iter = 0; iter < SFM_setting_window->ui.circle_iter_spinBox->value(); iter++)
		{
			oQProgressDialog.setWindowModality(Qt::ApplicationModal);
			oQProgressDialog.setWindowTitle(tr("Update iterations..."));
			oQProgressDialog.setRange(0, SFM_setting_window->ui.circle_iter_spinBox->value() + 1);
			oQProgressDialog.setMinimumDuration(0);
			oQProgressDialog.show();
			oQProgressDialog.setValue(iter + 1);
			ceres::Solver::Summary summary_ori;
			calculate_circle_thread thread_cal_circleori(xyz_sfm_ori, Keypoints_contours_in_queue, Keypoints_in_queue, Result_Ori,
				camQvec, camTvec, K, Dis_K, Dis_P, Dis_T, &summary_ori
				,  SFM_setting_window->ui.max_circle_iter_spinBox->value()
				, ((double)SFM_setting_window->ui.circle_stop_spinBox->value())* pow(10, SFM_setting_window->ui.circle_stop_spinBox_2->value())
				, SFM_setting_window->ui.thread_spinBox->value()
				, SFM_setting_window->ui.timeout_spinBox->value()
				, SFM_setting_window->ui.circle_loss_comboBox->currentIndex(), SFM_setting_window->ui.circle_loss_doubleSpinBox->value());
			thread_cal_circleori.start();
			while (1)
			{
				if (!thread_cal_circleori.isRunning())
				{
					break;
				}
				QString inf_log = "";
				for (int ss = summary_ori.iterations.size() - 1; ss >= 0; ss--)
				{
					inf_log = inf_log + QString("%1").arg(ss + 1, 4, 10, QLatin1Char('0')) + tr(". ");
					inf_log = inf_log + tr("[Cost]:") + QString::number(summary_ori.iterations[ss].cost, 'e', 3) + tr("; ");
					inf_log = inf_log + tr("[Gradient Norm]:") + QString::number(summary_ori.iterations[ss].gradient_norm, 'e', 3) + tr("; ");
					inf_log = inf_log + tr("[Step Norm]:") + QString::number(summary_ori.iterations[ss].step_norm, 'e', 3) + tr("; ");
					inf_log = inf_log + tr("[Iteration Time(second)]:") + QString::number(summary_ori.iterations[ss].iteration_time_in_seconds, 'e', 3) + tr(";");
					inf_log = inf_log + tr("\n");
				}
				cam_log_ui->ui.textEdit->setText(inf_log);
				loop_sleep(20);
			}
			Result_Ori = thread_cal_circleori.ori_serial;
			xyz_sfm_ori= thread_cal_circleori.point_clouds;

			std::vector<Eigen::Matrix<double, 3, 4>> K_serial;
			Eigen::Matrix<double, 3, 3> I_matrix = Eigen::Matrix<double, 3, 3>::Zero();
			I_matrix(0, 0) = K[0];
			I_matrix(0, 1) = K[4];
			I_matrix(0, 2) = K[2];
			I_matrix(1, 1) = K[1];
			I_matrix(1, 2) = K[3];
			I_matrix(2, 2) = 1;
			for (int ii = 0; ii < camQvec.size(); ii++)
			{
				Eigen::Matrix<double, 3, 3> R_matrix = camQvec[ii].matrix();
				Eigen::Matrix<double, 3, 4> RT_matrix = Eigen::Matrix<double, 3, 4>::Zero();
				RT_matrix(0, 0) = R_matrix(0, 0); RT_matrix(0, 1) = R_matrix(0, 1); RT_matrix(0, 2) = R_matrix(0, 2);
				RT_matrix(1, 0) = R_matrix(1, 0); RT_matrix(1, 1) = R_matrix(1, 1); RT_matrix(1, 2) = R_matrix(1, 2);
				RT_matrix(2, 0) = R_matrix(2, 0); RT_matrix(2, 1) = R_matrix(2, 1); RT_matrix(2, 2) = R_matrix(2, 2);
				RT_matrix(0, 3) = camTvec[ii].x(); RT_matrix(1, 3) = camTvec[ii].y(); RT_matrix(2, 3) = camTvec[ii].z();
				K_serial.push_back(I_matrix * RT_matrix);
			}
			QThreadPool::globalInstance()->setMaxThreadCount(SFM_setting_window->ui.thread_spinBox->value());

			std::vector<int> thread_for_circle_update_calcu_index;
			std::vector<bool> thread_for_circle_update_calcu_finished;
			std::vector<update_coded_circle_pixel_thread*> thread_for_circle_update_calcu;
			for (int ii = 0; ii < Keypoints_in_queue.size(); ii++)
			//for (int ii = 81; ii < 82; ii++)
			{
				thread_for_circle_update_calcu_index.push_back(ii);
				thread_for_circle_update_calcu_finished.push_back(false);
				thread_for_circle_update_calcu.push_back(new update_coded_circle_pixel_thread(K_serial[ii], camQvec[ii].matrix(), camTvec[ii], Keypoints_in_queue[ii],
					Result_Ori,xyz_sfm_ori, Keypoints_contours_in_queue[ii],
					Image_serial_name[Keypoints_index_in_queue[ii]], SFM_setting_window->ui.r_k_doubleSpinBox->value(),
					SFM_setting_window->ui.r_k_1_doubleSpinBox->value(), SFM_setting_window->ui.r_k_2_doubleSpinBox->value(),
					SFM_setting_window->ui.min_radius_doubleSpinBox->value(), SFM_setting_window->ui.max_radius_doubleSpinBox->value(),
					SFM_setting_window->ui.max_radius_error_doubleSpinBox->value(), mark_type, circle_code_type,
					contour_method, sub_pixel_method, SFM_setting_window->ui.max_aspect_ratio_doubleSpinBox->value()
					, SFM_setting_window->ui.min_contour_points_spinBox->value(), SFM_setting_window->ui.min_contour_number_spinBox->value()
					, SFM_setting_window->ui.min_gap_spinBox->value()
					, SFM_setting_window->ui.max_foreground_std_spinBox->value(), SFM_setting_window->ui.max_background_std_spinBox->value()));
				thread_for_circle_update_calcu[thread_for_circle_update_calcu.size() - 1]->setAutoDelete(false);
			}

			oQProgressDialog.setWindowModality(Qt::ApplicationModal);
			oQProgressDialog.setWindowTitle(tr("Update coded circle points..."));
			oQProgressDialog.setRange(0, thread_for_circle_update_calcu.size() + 1);
			oQProgressDialog.setMinimumDuration(0);
			oQProgressDialog.show();
			for (int ii = 0; ii < thread_for_circle_update_calcu.size(); ii++)
			{
				QThreadPool::globalInstance()->start(thread_for_circle_update_calcu[ii]);
			}
			oQProgressDialog.setValue(1);
			qApp->processEvents();

			while (1)
			{
				loop_sleep(10);
				int finished_number = 0;
				if (oQProgressDialog.wasCanceled())
				{
					KeyPoint_code_serial_update = KeyPoint_code_serial;
					QMessageBox::information(this, tr("Process"), tr("User Stop"), QMessageBox::Ok);
					in_calculation = false;
					return;
				}
				for (int ii = 0; ii < thread_for_circle_update_calcu.size(); ii++)
				{
					if (thread_for_circle_update_calcu_finished[ii])
					{
						finished_number++;
						continue;
					}
					if (thread_for_circle_update_calcu[ii]->is_finish)
					{
						finished_number++;
						Keypoints_in_queue[ii]= thread_for_circle_update_calcu[ii]->key_points_thread;
						KeyPoint_code_serial_update[Keypoints_index_in_queue[ii]]= thread_for_circle_update_calcu[ii]->key_points_thread;
						thread_for_circle_update_calcu_finished[ii] = true;
						cam_log_ui->ui.textEdit->setText(str_inf_for_log);
						delete thread_for_circle_update_calcu[ii];
					}
				}

				oQProgressDialog.setLabelText(QString::number(finished_number) + tr(" images have been update."));
				oQProgressDialog.setValue(finished_number + 1);
				qApp->processEvents();
				if (finished_number == thread_for_circle_update_calcu.size())
				{
					break;
				}
			}
			//for (int ii = 0; ii < 1; ii++)
			//{
			//	//std::cout << sss++ << std::endl;
			//	std::vector<Coded_detect_inf> Keypoints_in_serial_copy;
			//	QString img_name = Image_serial_name[Keypoints_index_in_queue[ii]];
			//	cv::Mat image_now = cv::imread(img_name.toLocal8Bit().toStdString());
			//	if (!image_now.data)
			//	{
			//		continue;
			//	}
			//	if (image_now.channels() == 3)
			//	{
			//		cv::cvtColor(image_now, image_now, CV_BGR2GRAY);
			//	}
			//	for (int jj = 0; jj < Keypoints_in_queue[ii].size(); jj++)
			//	{
			//		int find_sfm_flag = -1;
			//		for (int ss = 0; ss < xyz_sfm_ori.size(); ss++)
			//		{
			//			if (Keypoints_in_queue[ii][jj].code_num == xyz_sfm_ori[ss].code)
			//			{
			//				find_sfm_flag = ss;
			//				break;
			//			}
			//		}
			//		if (find_sfm_flag == -1)
			//		{
			//			continue;
			//		}
			//		double min_x = std::numeric_limits<double>::max();
			//		double min_y = std::numeric_limits<double>::max();
			//		double max_x = -std::numeric_limits<double>::max();
			//		double max_y = -std::numeric_limits<double>::max();
			//		for (int pp = 0; pp < Keypoints_contours_in_queue[ii][jj].size(); pp++)
			//		{
			//			min_x = min_x < Keypoints_contours_in_queue[ii][jj][pp].x ? min_x : Keypoints_contours_in_queue[ii][jj][pp].x;
			//			min_y = min_y < Keypoints_contours_in_queue[ii][jj][pp].y ? min_y : Keypoints_contours_in_queue[ii][jj][pp].y;
			//			max_x = max_x > Keypoints_contours_in_queue[ii][jj][pp].x ? max_x : Keypoints_contours_in_queue[ii][jj][pp].x;
			//			max_y = max_y > Keypoints_contours_in_queue[ii][jj][pp].y ? max_y : Keypoints_contours_in_queue[ii][jj][pp].y;
			//		}
			//		double muli_cof = 1.25;
			//		min_x = Keypoints_in_queue[ii][jj].x - muli_cof * SFM_setting_window->ui.r_k_2_doubleSpinBox->value() * (Keypoints_in_queue[ii][jj].x - min_x);
			//		max_x = Keypoints_in_queue[ii][jj].x + muli_cof * SFM_setting_window->ui.r_k_2_doubleSpinBox->value() * (-Keypoints_in_queue[ii][jj].x + max_x);
			//		min_y = Keypoints_in_queue[ii][jj].y - muli_cof * SFM_setting_window->ui.r_k_2_doubleSpinBox->value() * (Keypoints_in_queue[ii][jj].y - min_y);
			//		max_y = Keypoints_in_queue[ii][jj].y + muli_cof * SFM_setting_window->ui.r_k_2_doubleSpinBox->value() * (-Keypoints_in_queue[ii][jj].y + max_y);
			//		min_x = floor(min_x);
			//		min_y = floor(min_y);
			//		max_x = ceil(max_x);
			//		max_y = ceil(max_y);
			//		min_x = min_x < 0 ? 0 : min_x;
			//		min_y = min_y < 0 ? 0 : min_y;
			//		max_x = max_x >= image_now.cols ? (image_now.cols - 1) : max_x;
			//		max_y = max_y >= image_now.rows ? (image_now.rows - 1) : max_y;
			//		cv::Mat Circle_image_Old = image_now(cv::Range(min_y, max_y), cv::Range(min_x, max_x));
			//		int window_R_new = ceil(muli_cof * ((max_y - min_y) > (max_x - min_x) ? (max_y - min_y) : (max_x - min_x)));
			//		//Eigen::MatrixXd new_ori_temp = camQvec[ii].toRotationMatrix() * Result_Ori[ii].matrix();
			//		Eigen::MatrixXd new_ori_temp = Result_Ori[find_sfm_flag].matrix();
			//		Eigen::Vector3d new_ori(new_ori_temp(0, 0), new_ori_temp(1, 0), new_ori_temp(2, 0));
			//		Eigen::Vector3d cros_i(1, 0, 0);
			//		Eigen::Vector3d cros_j(0, 1, 0);
			//		Eigen::Vector3d cros_a, cros_b;
			//		cros_a = new_ori.cross(cros_i);
			//		if (cros_a.x() == 0 && cros_a.y() == 0 && cros_a.z() == 0)
			//		{
			//			cros_a = new_ori.cross(cros_j);
			//		}
			//		cros_b = new_ori.cross(cros_a);
			//		cros_a.normalize();
			//		cros_b.normalize();
			//		Eigen::MatrixXd new_sfm_xyz = Eigen::Vector3d(xyz_sfm_ori[find_sfm_flag].x, xyz_sfm_ori[find_sfm_flag].y, xyz_sfm_ori[find_sfm_flag].z).matrix();
			//		Eigen::Vector4d old_point_1_4d(new_sfm_xyz(0, 0) + muli_cof * SFM_setting_window->ui.r_k_2_doubleSpinBox->value() * xyz_sfm_ori[find_sfm_flag].R * cros_a(0),
			//			new_sfm_xyz(1, 0) + muli_cof * SFM_setting_window->ui.r_k_2_doubleSpinBox->value() * xyz_sfm_ori[find_sfm_flag].R * cros_a(1), new_sfm_xyz(2, 0) + muli_cof * SFM_setting_window->ui.r_k_2_doubleSpinBox->value() * xyz_sfm_ori[find_sfm_flag].R * cros_a(2), 1);
			//		Eigen::Vector4d old_point_2_4d(new_sfm_xyz(0, 0) + muli_cof * SFM_setting_window->ui.r_k_2_doubleSpinBox->value() * xyz_sfm_ori[find_sfm_flag].R * cros_b(0),
			//			new_sfm_xyz(1, 0) + muli_cof * SFM_setting_window->ui.r_k_2_doubleSpinBox->value() * xyz_sfm_ori[find_sfm_flag].R * cros_b(1), new_sfm_xyz(2, 0) + muli_cof * SFM_setting_window->ui.r_k_2_doubleSpinBox->value() * xyz_sfm_ori[find_sfm_flag].R * cros_b(2), 1);
			//		Eigen::Vector4d old_point_3_4d(new_sfm_xyz(0, 0) - muli_cof * SFM_setting_window->ui.r_k_2_doubleSpinBox->value() * xyz_sfm_ori[find_sfm_flag].R * cros_a(0),
			//			new_sfm_xyz(1, 0) - muli_cof * SFM_setting_window->ui.r_k_2_doubleSpinBox->value() * xyz_sfm_ori[find_sfm_flag].R * cros_a(1), new_sfm_xyz(2, 0) - muli_cof * SFM_setting_window->ui.r_k_2_doubleSpinBox->value() * xyz_sfm_ori[find_sfm_flag].R * cros_a(2), 1);
			//		Eigen::Vector4d old_point_4_4d(new_sfm_xyz(0, 0) - muli_cof * SFM_setting_window->ui.r_k_2_doubleSpinBox->value() * xyz_sfm_ori[find_sfm_flag].R * cros_b(0),
			//			new_sfm_xyz(1, 0) - muli_cof * SFM_setting_window->ui.r_k_2_doubleSpinBox->value() * xyz_sfm_ori[find_sfm_flag].R * cros_b(1), new_sfm_xyz(2, 0) - muli_cof * SFM_setting_window->ui.r_k_2_doubleSpinBox->value() * xyz_sfm_ori[find_sfm_flag].R * cros_b(2), 1);
			//		Eigen::Vector3d old_point_1_3d = K_serial[ii] * old_point_1_4d.matrix();
			//		Eigen::Vector3d old_point_2_3d = K_serial[ii] * old_point_2_4d.matrix();
			//		Eigen::Vector3d old_point_3_3d = K_serial[ii] * old_point_3_4d.matrix();
			//		Eigen::Vector3d old_point_4_3d = K_serial[ii] * old_point_4_4d.matrix();
			//		cv::Point2f old_point_1(old_point_1_3d.x() / old_point_1_3d.z() - min_x, old_point_1_3d.y() / old_point_1_3d.z() - min_y);
			//		cv::Point2f old_point_2(old_point_2_3d.x() / old_point_2_3d.z() - min_x, old_point_2_3d.y() / old_point_2_3d.z() - min_y);
			//		cv::Point2f old_point_3(old_point_3_3d.x() / old_point_3_3d.z() - min_x, old_point_3_3d.y() / old_point_3_3d.z() - min_y);
			//		cv::Point2f old_point_4(old_point_4_3d.x() / old_point_4_3d.z() - min_x, old_point_4_3d.y() / old_point_4_3d.z() - min_y);
			//		if (jj == 1)
			//		{
			//			std::cout << old_point_1 << std::endl;
			//			std::cout << old_point_2 << std::endl;
			//			std::cout << old_point_3 << std::endl;
			//			std::cout << old_point_4 << std::endl;
			//		}
			//		cv::Point2f new_point_1(window_R_new, 0);
			//		cv::Point2f new_point_2(0, window_R_new);
			//		cv::Point2f new_point_3(window_R_new, 2 * window_R_new);
			//		cv::Point2f new_point_4(2 * window_R_new, window_R_new);
			//		cv::Point2f AffinePointsSrc[4] = { old_point_1, old_point_2, old_point_3, old_point_4 };
			//		cv::Point2f AffinePointsDst[4] = { new_point_1, new_point_2, new_point_3, new_point_4 };
			//		cv::Mat Circle_image_New;
			//		cv::Mat per_tf = cv::getPerspectiveTransform(AffinePointsSrc, AffinePointsDst);
			//		cv::warpPerspective(Circle_image_Old, Circle_image_New, per_tf, cv::Size(2 * window_R_new + 1, 2 * window_R_new + 1), CV_INTER_CUBIC);
			//		cv::Mat code_point_temp;
			//		std::vector<std::vector<cv::Point>> pixels;
			//		//cv::namedWindow("cannys", cv::WINDOW_NORMAL);
			//		//cv::imshow("cannys", Circle_image_Old);
			//		//cv::waitKey(0);
			//		//cv::namedWindow("cannyqqq", cv::WINDOW_NORMAL);
			//		//cv::imshow("cannyqqq", Circle_image_New);
			//		//cv::waitKey(0);
			//		ImageDetectMethod::detectcodecircle(Circle_image_New, code_point_temp, pixels, SFM_setting_window->ui.r_k_doubleSpinBox->value(),
			//			SFM_setting_window->ui.r_k_1_doubleSpinBox->value(), SFM_setting_window->ui.r_k_2_doubleSpinBox->value(),
			//			2, window_R_new * 2 + 1,
			//			SFM_setting_window->ui.max_radius_error_doubleSpinBox->value(), mark_type, circle_code_type,
			//			contour_method, sub_pixel_method, SFM_setting_window->ui.max_aspect_ratio_doubleSpinBox->value()
			//			, SFM_setting_window->ui.min_contour_points_spinBox->value(), SFM_setting_window->ui.min_contour_number_spinBox->value()
			//			, SFM_setting_window->ui.min_gap_spinBox->value()
			//			, SFM_setting_window->ui.max_foreground_std_spinBox->value(), SFM_setting_window->ui.max_background_std_spinBox->value());
			//		bool can_be_ex = false;
			//		for (unsigned int gg = 0; gg < code_point_temp.rows; gg++)
			//		{
			//			if (code_point_temp.at<float>(gg, 0) == Keypoints_in_queue[ii][jj].code_num)
			//			{
			//				can_be_ex = true;
			//				cv::Mat_<double> mat_pt(3, 1);
			//				mat_pt(0, 0) = code_point_temp.at<float>(gg, 1);
			//				mat_pt(1, 0) = code_point_temp.at<float>(gg, 2);
			//				mat_pt(2, 0) = 1;
			//				cv::Mat mat_pt_view = per_tf.inv() * mat_pt;
			//				double a1 = mat_pt_view.at<double>(0, 0);
			//				double a2 = mat_pt_view.at<double>(1, 0);
			//				double a3 = mat_pt_view.at<double>(2, 0);
			//				double new_x = a1 * 1.0 / a3 + min_x;
			//				double new_y = a2 * 1.0 / a3 + min_y;
			//				Keypoints_in_queue[ii][jj].x = new_x;
			//				Keypoints_in_queue[ii][jj].y = new_y;
			//				KeyPoint_code_serial[Keypoints_index_in_queue[ii]][jj].x = new_x;
			//				KeyPoint_code_serial[Keypoints_index_in_queue[ii]][jj].y = new_y;
			//				Keypoints_in_serial_copy.push_back(Keypoints_in_queue[ii][jj]);
			//				break;
			//			}
			//		}
			//	}
			//}
			if (SFM_setting_window->ui.optimize_checkBox->isChecked())
			{
				K[0] = fx;
				K[1] = fy;
				K[2] = cx;
				K[3] = cy;
				K[4] = 0;
				Dis_K[0] = 0;
				Dis_K[1] = 0;
				Dis_K[2] = 0;
				Dis_K[3] = 0;
				Dis_K[4] = 0;
				Dis_K[5] = 0;
				Dis_P[0] = 0;
				Dis_P[1] = 0;
				Dis_T[0] = 0;
				Dis_T[1] = 0;
				Dis_T[2] = 0;
				Dis_T[3] = 0;
			}
			need_optimize_C = true;
			if (SFM_setting_window->ui.fixed_c_checkBox->isChecked())
			{
				need_optimize_C = false;
			}
			need_optimize_F = true;
			if (SFM_setting_window->ui.fixed_f_checkBox->isChecked())
			{
				need_optimize_F = false;
			}
			ceres::Solver::Summary summary_sfm;
			optimize_SFM_incremental_thread thread_cal_new(Keypoints_in_queue, xyz_sfm_ori, camQvec, camTvec, Re_project_Map, K, Dis_K, Dis_P, Dis_T, &summary_sfm
				, SFM_setting_window->ui.max_iter_spinBox->value()
				, ((double)SFM_setting_window->ui.stop_spinBox->value()) * pow(10, SFM_setting_window->ui.stop_spinBox_2->value())
				, SFM_setting_window->ui.thread_spinBox->value()
				, SFM_setting_window->ui.timeout_spinBox->value()
				, SFM_setting_window->ui.loss_comboBox->currentIndex(), SFM_setting_window->ui.loss_doubleSpinBox->value()
				, SFM_setting_window->ui.k_distort_comboBox->currentText().toInt()
				, SFM_setting_window->ui.p_distort_comboBox->currentText().toInt()
				, SFM_setting_window->ui.s_distort_comboBox->currentText().toInt()
				, SFM_setting_window->ui.fixed_camera_checkBox->isChecked()
				, need_optimize_F
				, need_optimize_C
				, SFM_setting_window->ui.skew_checkBox->isChecked()
				, SFM_setting_window->ui.uniform_f_checkBox->isChecked(), range_f, range_c);
			thread_cal_new.start();
			while (1)
			{
				if (!thread_cal_new.isRunning())
				{
					break;
				}
				QString inf_log = "";
				for (int ss = summary_sfm.iterations.size() - 1; ss >= 0; ss--)
				{
					inf_log = inf_log + QString("%1").arg(ss + 1, 4, 10, QLatin1Char('0')) + tr(". ");
					inf_log = inf_log + tr("[Cost]:") + QString::number(summary_sfm.iterations[ss].cost, 'e', 3) + tr("; ");
					inf_log = inf_log + tr("[Gradient Norm]:") + QString::number(summary_sfm.iterations[ss].gradient_norm, 'e', 3) + tr("; ");
					inf_log = inf_log + tr("[Step Norm]:") + QString::number(summary_sfm.iterations[ss].step_norm, 'e', 3) + tr("; ");
					inf_log = inf_log + tr("[Iteration Time(second)]:") + QString::number(summary_sfm.iterations[ss].iteration_time_in_seconds, 'e', 3) + tr(";");
					inf_log = inf_log + tr("\n");
				}
				cam_log_ui->ui.textEdit->setText(inf_log);
				loop_sleep(20);
			}

			xyz_sfm_ori = thread_cal_new.point_clouds_thread;
			camQvec = thread_cal_new.camQvec_thread;
			camTvec = thread_cal_new.camTvec_thread;
			Re_project_Map = thread_cal_new.Re_project_Map_thread;
			K = thread_cal_new.camK_thread;
			Dis_K = thread_cal_new.camDK_thread;
			Dis_P = thread_cal_new.camDP_thread;
			Dis_T = thread_cal_new.camDT_thread;
			str_inf_for_log =
				tr("Global BA completion.\n") + str_inf_for_log;
			cam_log_ui->ui.textEdit->setText(str_inf_for_log);
			Result_Camera_K[0] = K[0];
			Result_Camera_K[1] = K[1];
			Result_Camera_K[2] = K[2];
			Result_Camera_K[3] = K[3];
			Result_Camera_K[4] = K[4];
			Result_Camera_Dis_K[0] = Dis_K[0];
			Result_Camera_Dis_K[1] = Dis_K[1];
			Result_Camera_Dis_K[2] = Dis_K[2];
			Result_Camera_Dis_K[3] = Dis_K[3];
			Result_Camera_Dis_K[4] = Dis_K[4];
			Result_Camera_Dis_K[5] = Dis_K[5];
			Result_Camera_Dis_P[0] = Dis_P[0];
			Result_Camera_Dis_P[1] = Dis_P[1];
			Result_Camera_Dis_T[0] = Dis_T[0];
			Result_Camera_Dis_T[1] = Dis_T[1];
			Result_Camera_Dis_T[2] = Dis_T[2];
			Result_Camera_Dis_T[3] = Dis_T[3];
			for (int pp = 1; pp < Keypoints_R_in_queue.size(); pp++)
			{
				cv::Mat mat_temp;
				cv::eigen2cv(camQvec[pp].matrix(), mat_temp);
				Keypoints_R_in_queue[pp] = mat_temp;
				cv::Mat mat_temp2;
				cv::eigen2cv(camTvec[pp], mat_temp2);
				Keypoints_T_in_queue[pp] = mat_temp2;
			}
			Result_err_every.clear();
			Result_R.clear();
			Result_T.clear();
			for (int ii = 0; ii < KeyPoint_code_serial.size(); ii++)
			{
				bool exist = false;
				for (int jj = 0; jj < Keypoints_index_in_queue.size(); jj++)
				{
					if (ii == Keypoints_index_in_queue[jj])
					{
						exist = true;
						Result_err_every.push_back(Re_project_Map[jj]);
						Result_R.push_back(camQvec[jj]);
						Result_T.push_back(camTvec[jj]);
						break;
					}
				}
				if (!exist)
				{
					std::vector<cv::Point2d> temp_err = {};
					for (int jj = 0; jj < KeyPoint_code_serial[ii].size(); jj++)
					{
						temp_err.push_back(cv::Point2d(std::numeric_limits<double>::max(), std::numeric_limits<double>::max()));
					}
					Result_err_every.push_back(temp_err);
					Result_R.push_back(Eigen::Quaterniond(0, 0, 0, 0));
					Result_T.push_back(Eigen::Vector3d(0, 0, 0));
				}
			}

			K_camera = cv::Mat(cv::Matx33d(
				K[0], K[4], K[2],
				0, K[1], K[3],
				0, 0, 1));
			dis_coff.at<double>(0) = Dis_K[0];
			dis_coff.at<double>(1) = Dis_K[1];
			dis_coff.at<double>(2) = Dis_P[0];
			dis_coff.at<double>(3) = Dis_P[1];
			dis_coff.at<double>(4) = Dis_K[2];
			dis_coff.at<double>(5) = Dis_K[3];
			dis_coff.at<double>(6) = Dis_K[4];
			dis_coff.at<double>(7) = Dis_K[5];
			dis_coff.at<double>(8) = Dis_T[0];
			dis_coff.at<double>(9) = Dis_T[2];
			dis_coff.at<double>(10) = Dis_T[1];
			dis_coff.at<double>(11) = Dis_T[3];
		}
	}

	ui.result_fx_doubleSpinBox->setValue(K[0]);
	ui.result_fy_doubleSpinBox->setValue(K[1]);
	ui.result_cx_doubleSpinBox->setValue(K[2]);
	ui.result_cy_doubleSpinBox->setValue(K[3]);
	ui.result_shear_doubleSpinBox->setValue(K[4]);
	ui.result_dis_K_doubleSpinBox_1->setValue(Dis_K[0]);
	ui.result_dis_K_doubleSpinBox_2->setValue(Dis_K[1]);
	ui.result_dis_K_doubleSpinBox_3->setValue(Dis_K[2]);
	ui.result_dis_K_doubleSpinBox_4->setValue(Dis_K[3]);
	ui.result_dis_K_doubleSpinBox_5->setValue(Dis_K[4]);
	ui.result_dis_K_doubleSpinBox_6->setValue(Dis_K[5]);
	ui.result_dis_P_doubleSpinBox->setValue(Dis_P[0]);
	ui.result_dis_P_doubleSpinBox_2->setValue(Dis_P[1]);
	ui.result_dis_T_doubleSpinBox_1->setValue(Dis_T[0]);
	ui.result_dis_T_doubleSpinBox_2->setValue(Dis_T[1]);
	ui.result_dis_T_doubleSpinBox_3->setValue(Dis_T[2]);
	ui.result_dis_T_doubleSpinBox_4->setValue(Dis_T[3]);
	if (range_f != nullptr)
	{
		delete[] range_f;
		range_f = nullptr;
	}	
	if (range_c != nullptr)
	{
		delete[] range_c;
		range_c = nullptr;
	}
	cam_log_ui->setWindowOpacity(1);
	oQProgressDialog.cancel();
	Result_xyz_serial = xyz_sfm_ori;
	xyz_need_levelling_serial.clear();
	xyz_need_levelling_code_serial.clear();
	for (int qs = 0; qs < Result_xyz_serial.size(); qs++)
	{
		xyz_need_levelling_serial.push_back(false);
		xyz_need_levelling_code_serial.push_back(Result_xyz_serial[qs].code);
	}
	std::sort(xyz_need_levelling_code_serial.begin(), xyz_need_levelling_code_serial.end());

	update_show_view();
	need_recal_keypoints = false;
	in_calculation = false;

	update_points_view(); 
	update_table_view();
}
double SFM_module::kahanSum_double(std::vector<double> nums)
{
	double sum = 0.0;
	double c = 0.0;
	for (auto num : nums) {
		double y = num - c;
		double t = sum + y;
		c = (t - sum) - y;
		sum = t;
	}
	return sum;
}
float SFM_module::kahanSum_float(std::vector<float> nums)
{
	float sum = 0.0f;
	float c = 0.0f;
	for (auto num : nums) {
		float y = num - c;
		float t = sum + y;
		c = (t - sum) - y;
		sum = t;
	}
	return sum;
}
ellipse_view_pose SFM_module::calculate_ellipse(std::vector<cv::Point> ellipse_contours,
	cv::Mat K, cv::Mat Dis, Eigen::Quaterniond R, Eigen::Vector3d T)
{
	std::vector<cv::Point2d> ellipse_contours_copy;
	for (int ii = 0; ii < ellipse_contours.size(); ii++)
	{
		ellipse_contours_copy.push_back(cv::Point2d(ellipse_contours[ii].x, ellipse_contours[ii].y));
	}
	cv::undistortPoints(ellipse_contours_copy, ellipse_contours_copy, K, Dis, K);

	double Ell_A, Ell_B, Ell_C, Ell_D, Ell_E, Ell_F;
	double sum_x = 0, sum_y = 0, sum_xx = 0, sum_xy = 0, sum_yy = 0, sum_xxx = 0, sum_xxy = 0, sum_xyy = 0, sum_yyy = 0, sum_xxxy = 0, sum_xxyy = 0, sum_xyyy = 0, sum_yyyy = 0;
	std::vector<double> c_x, c_y, c_xx, c_xy, c_yy, c_xxx, c_xxy, c_xyy, c_yyy, c_xxxy, c_xxyy, c_xyyy, c_yyyy;
	double max_x = -std::numeric_limits<double>::max(), max_y = -std::numeric_limits<double>::max();
	for (int ii = 0; ii < ellipse_contours_copy.size(); ii++)
	{
		max_x = max_x < abs(ellipse_contours_copy[ii].x) ? abs(ellipse_contours_copy[ii].x) : max_x;
		max_y = max_y < abs(ellipse_contours_copy[ii].y) ? abs(ellipse_contours_copy[ii].y) : max_y;
	}
	for (int ii = 0; ii < ellipse_contours_copy.size(); ii++)
	{
		ellipse_contours_copy[ii].x /= max_x;
		ellipse_contours_copy[ii].y /= max_y;
		c_x.push_back(ellipse_contours_copy[ii].x);
		c_y.push_back(ellipse_contours_copy[ii].y);
		c_xx.push_back(ellipse_contours_copy[ii].x * ellipse_contours_copy[ii].x);
		c_xy.push_back(ellipse_contours_copy[ii].x * ellipse_contours_copy[ii].y);
		c_yy.push_back(ellipse_contours_copy[ii].y * ellipse_contours_copy[ii].y);
		c_xxx.push_back(ellipse_contours_copy[ii].x * ellipse_contours_copy[ii].x * ellipse_contours_copy[ii].x);
		c_xxy.push_back(ellipse_contours_copy[ii].x * ellipse_contours_copy[ii].x * ellipse_contours_copy[ii].y);
		c_xyy.push_back(ellipse_contours_copy[ii].x * ellipse_contours_copy[ii].y * ellipse_contours_copy[ii].y);
		c_yyy.push_back(ellipse_contours_copy[ii].y * ellipse_contours_copy[ii].y * ellipse_contours_copy[ii].y);
		c_xxxy.push_back(ellipse_contours_copy[ii].x * ellipse_contours_copy[ii].x * ellipse_contours_copy[ii].x * ellipse_contours_copy[ii].y);
		c_xxyy.push_back(ellipse_contours_copy[ii].x * ellipse_contours_copy[ii].x * ellipse_contours_copy[ii].y * ellipse_contours_copy[ii].y);
		c_xyyy.push_back(ellipse_contours_copy[ii].x * ellipse_contours_copy[ii].y * ellipse_contours_copy[ii].y * ellipse_contours_copy[ii].y);
		c_yyyy.push_back(ellipse_contours_copy[ii].y * ellipse_contours_copy[ii].y * ellipse_contours_copy[ii].y * ellipse_contours_copy[ii].y);
	}
	sum_x = kahanSum_double(c_x);
	sum_y = kahanSum_double(c_y);
	sum_xx = kahanSum_double(c_xx);
	sum_xy = kahanSum_double(c_xy);
	sum_yy = kahanSum_double(c_yy);
	sum_xxx = kahanSum_double(c_xxx);
	sum_xxy = kahanSum_double(c_xxy);
	sum_xyy = kahanSum_double(c_xyy);
	sum_yyy = kahanSum_double(c_yyy);
	sum_xxxy = kahanSum_double(c_xxxy);
	sum_xxyy = kahanSum_double(c_xxyy);
	sum_xyyy = kahanSum_double(c_xyyy);
	sum_yyyy = kahanSum_double(c_yyyy);
	Eigen::MatrixXd B = Eigen::MatrixXd::Zero(5, 5);
	B(0, 0) = sum_xxyy; B(0, 1) = sum_xyyy; B(0, 2) = sum_xxy; B(0, 3) = sum_xyy; B(0, 4) = sum_xy;
	B(1, 0) = sum_xyyy; B(1, 1) = sum_yyyy; B(1, 2) = sum_xyy; B(1, 3) = sum_yyy; B(1, 4) = sum_yy;
	B(2, 0) = sum_xxy; B(2, 1) = sum_xyy; B(2, 2) = sum_xx; B(2, 3) = sum_xy; B(2, 4) = sum_x;
	B(3, 0) = sum_xyy; B(3, 1) = sum_yyy; B(3, 2) = sum_xy; B(3, 3) = sum_yy; B(3, 4) = sum_y;
	B(4, 0) = sum_xy; B(4, 1) = sum_yy; B(4, 2) = sum_x; B(4, 3) = sum_y; B(4, 4) = ellipse_contours_copy.size();
	Eigen::MatrixXd C = Eigen::MatrixXd::Zero(5, 1);
	C(0, 0) = sum_xxxy; C(1, 0) = sum_xxyy; C(2, 0) = sum_xxx; C(3, 0) = sum_xxy; C(4, 0) = sum_xx;
	Eigen::MatrixXd D = -B.inverse() * C;
	Ell_A = 1.0 / (max_x * max_x * D(4));
	Ell_B = D(0) / (max_x * max_y * D(4));
	Ell_C = D(1) / (max_y * max_y * D(4));
	Ell_D = D(2) / (max_x * D(4));
	Ell_E = D(3) / (max_y * D(4));
	Ell_F = 1.0;
	//cv::RotatedRect rRect = cv::fitEllipseAMS(ellipse_contours_copy);
	//rRect.angle = rRect.angle* M_PI / 180;
	//double Ell_A = (rRect.size.height * rRect.size.height / 4.0 * cos(rRect.angle) * cos(rRect.angle) +
	//		rRect.size.width * rRect.size.width / 4.0 * sin(rRect.angle) * sin(rRect.angle));
	//double Ell_B = - 2 * (rRect.size.width * rRect.size.width - rRect.size.height * rRect.size.height) / 4.0 * sin(rRect.angle) * cos(rRect.angle);
	//double Ell_C = (rRect.size.height * rRect.size.height / 4.0 * sin(rRect.angle) * sin(rRect.angle) +
	//	rRect.size.width * rRect.size.width / 4.0 * cos(rRect.angle) * cos(rRect.angle));
	//double Ell_D = - 2 * rRect.center.x * Ell_A + rRect.center.y * Ell_B;
	//double Ell_E = - 2 * rRect.center.y * Ell_C + rRect.center.x * Ell_B;
	//double Ell_F = Ell_A * rRect.center.x * rRect.center.x + rRect.center.y * rRect.center.y * Ell_B + rRect.center.x * rRect.center.y * Ell_C -
	//	-rRect.size.width * rRect.size.width * rRect.size.height * rRect.size.height / 16.0;

	//Ell_A /= Ell_F;
	//Ell_B /= Ell_F;
	//Ell_C /= Ell_F;
	//Ell_D /= Ell_F;
	//Ell_E /= Ell_F;
	//Ell_F = 1.0;

	Eigen::MatrixXd Ell_g_matrix(3, 3);
	Ell_g_matrix(0, 0) = Ell_A;
	Ell_g_matrix(0, 1) = Ell_B / 2.0;
	Ell_g_matrix(0, 2) = Ell_D / 2.0;
	Ell_g_matrix(1, 0) = Ell_B / 2.0;
	Ell_g_matrix(1, 1) = Ell_C;
	Ell_g_matrix(1, 2) = Ell_E / 2.0;
	Ell_g_matrix(2, 0) = Ell_D / 2.0;
	Ell_g_matrix(2, 1) = Ell_E / 2.0;
	Ell_g_matrix(2, 2) = Ell_F;


	Eigen::MatrixXd Ins_matrix(3, 3);
	cv::cv2eigen(K, Ins_matrix);
	Eigen::MatrixXd Rt_matrix= Eigen::MatrixXd::Zero(3,4);
	Rt_matrix(0, 0) = 1;
	Rt_matrix(1, 1) = 1;
	Rt_matrix(2, 2) = 1;
	//for (int ii = 0; ii < 3; ii++)
	//{
	//	for (int jj = 0; jj < 3; jj++)
	//	{
	//		Rt_matrix(ii, jj) = R.matrix()(ii, jj);
	//	}
	//	Rt_matrix(ii, 3) = T(ii);
	//}
	Eigen::MatrixXd K_matrix = Ins_matrix * Rt_matrix;
	Eigen::MatrixXd Ell_q_matrix= K_matrix.transpose() * Ell_g_matrix * K_matrix;
	Eigen::MatrixXd Ell_Q_matrix = Ell_q_matrix.block<3, 3>(0, 0);

	Eigen::EigenSolver<Eigen::MatrixXd> es(Ell_Q_matrix);


	double eigenvalue_1_temp = es.eigenvalues()(0).real();
	double eigenvalue_2_temp = es.eigenvalues()(1).real();
	double eigenvalue_3_temp = es.eigenvalues()(2).real();
	double eigenvalue_1, eigenvalue_2, eigenvalue_3;
	Eigen::Vector3d eigenvector_1;
	Eigen::Vector3d eigenvector_2;
	Eigen::Vector3d eigenvector_3;

	if (eigenvalue_1_temp * eigenvalue_2_temp >= 0)
	{
		if (abs(eigenvalue_1_temp) > abs(eigenvalue_2_temp))
		{
			eigenvalue_1 = eigenvalue_1_temp;
			eigenvalue_2 = eigenvalue_2_temp;
			eigenvalue_3 = eigenvalue_3_temp;
			int index_1 = 0, index_2 = 1, index_3 = 2;
			eigenvector_1(0) = es.eigenvectors()(0, index_1).real();
			eigenvector_1(1) = es.eigenvectors()(1, index_1).real();
			eigenvector_1(2) = es.eigenvectors()(2, index_1).real();
			eigenvector_2(0) = es.eigenvectors()(0, index_2).real();
			eigenvector_2(1) = es.eigenvectors()(1, index_2).real();
			eigenvector_2(2) = es.eigenvectors()(2, index_2).real();
			eigenvector_3(0) = es.eigenvectors()(0, index_3).real();
			eigenvector_3(1) = es.eigenvectors()(1, index_3).real();
			eigenvector_3(2) = es.eigenvectors()(2, index_3).real();
		}
		else
		{
			eigenvalue_1 = eigenvalue_2_temp;
			eigenvalue_2 = eigenvalue_1_temp;
			eigenvalue_3 = eigenvalue_3_temp;
			int index_1 = 1, index_2 = 0, index_3 = 2;
			eigenvector_1(0) = es.eigenvectors()(0, index_1).real();
			eigenvector_1(1) = es.eigenvectors()(1, index_1).real();
			eigenvector_1(2) = es.eigenvectors()(2, index_1).real();
			eigenvector_2(0) = es.eigenvectors()(0, index_2).real();
			eigenvector_2(1) = es.eigenvectors()(1, index_2).real();
			eigenvector_2(2) = es.eigenvectors()(2, index_2).real();
			eigenvector_3(0) = es.eigenvectors()(0, index_3).real();
			eigenvector_3(1) = es.eigenvectors()(1, index_3).real();
			eigenvector_3(2) = es.eigenvectors()(2, index_3).real();
		}
	}
	else if (eigenvalue_1_temp * eigenvalue_3_temp > 0)
	{
		if (abs(eigenvalue_1_temp) > abs(eigenvalue_3_temp))
		{
			eigenvalue_1 = eigenvalue_1_temp;
			eigenvalue_2 = eigenvalue_3_temp;
			eigenvalue_3 = eigenvalue_2_temp;
			int index_1 = 0, index_2 = 2, index_3 = 1;
			eigenvector_1(0) = es.eigenvectors()(0, index_1).real();
			eigenvector_1(1) = es.eigenvectors()(1, index_1).real();
			eigenvector_1(2) = es.eigenvectors()(2, index_1).real();
			eigenvector_2(0) = es.eigenvectors()(0, index_2).real();
			eigenvector_2(1) = es.eigenvectors()(1, index_2).real();
			eigenvector_2(2) = es.eigenvectors()(2, index_2).real();
			eigenvector_3(0) = es.eigenvectors()(0, index_3).real();
			eigenvector_3(1) = es.eigenvectors()(1, index_3).real();
			eigenvector_3(2) = es.eigenvectors()(2, index_3).real();
		}
		else
		{
			eigenvalue_1 = eigenvalue_3_temp;
			eigenvalue_2 = eigenvalue_1_temp;
			eigenvalue_3 = eigenvalue_2_temp;
			int index_1 = 2, index_2 = 0, index_3 = 1;
			eigenvector_1(0) = es.eigenvectors()(0, index_1).real();
			eigenvector_1(1) = es.eigenvectors()(1, index_1).real();
			eigenvector_1(2) = es.eigenvectors()(2, index_1).real();
			eigenvector_2(0) = es.eigenvectors()(0, index_2).real();
			eigenvector_2(1) = es.eigenvectors()(1, index_2).real();
			eigenvector_2(2) = es.eigenvectors()(2, index_2).real();
			eigenvector_3(0) = es.eigenvectors()(0, index_3).real();
			eigenvector_3(1) = es.eigenvectors()(1, index_3).real();
			eigenvector_3(2) = es.eigenvectors()(2, index_3).real();
		}
	}
	else
	{
		if (abs(eigenvalue_2_temp) > abs(eigenvalue_3_temp))
		{
			eigenvalue_1 = eigenvalue_2_temp;
			eigenvalue_2 = eigenvalue_3_temp;
			eigenvalue_3 = eigenvalue_1_temp;
			int index_1 = 1, index_2 = 2, index_3 = 0;
			eigenvector_1(0) = es.eigenvectors()(0, index_1).real();
			eigenvector_1(1) = es.eigenvectors()(1, index_1).real();
			eigenvector_1(2) = es.eigenvectors()(2, index_1).real();
			eigenvector_2(0) = es.eigenvectors()(0, index_2).real();
			eigenvector_2(1) = es.eigenvectors()(1, index_2).real();
			eigenvector_2(2) = es.eigenvectors()(2, index_2).real();
			eigenvector_3(0) = es.eigenvectors()(0, index_3).real();
			eigenvector_3(1) = es.eigenvectors()(1, index_3).real();
			eigenvector_3(2) = es.eigenvectors()(2, index_3).real();
		}
		else
		{
			eigenvalue_1 = eigenvalue_3_temp;
			eigenvalue_2 = eigenvalue_2_temp;
			eigenvalue_3 = eigenvalue_1_temp;
			int index_1 = 2, index_2 = 1, index_3 = 0;
			eigenvector_1(0) = es.eigenvectors()(0, index_1).real();
			eigenvector_1(1) = es.eigenvectors()(1, index_1).real();
			eigenvector_1(2) = es.eigenvectors()(2, index_1).real();
			eigenvector_2(0) = es.eigenvectors()(0, index_2).real();
			eigenvector_2(1) = es.eigenvectors()(1, index_2).real();
			eigenvector_2(2) = es.eigenvectors()(2, index_2).real();
			eigenvector_3(0) = es.eigenvectors()(0, index_3).real();
			eigenvector_3(1) = es.eigenvectors()(1, index_3).real();
			eigenvector_3(2) = es.eigenvectors()(2, index_3).real();
		}
	}

	if (eigenvalue_1_temp >=0 && eigenvalue_2_temp >= 0 && eigenvalue_3_temp >= 0)
	{
		double val_1 = abs(eigenvalue_1_temp - eigenvalue_2_temp);
		double val_2 = abs(eigenvalue_1_temp - eigenvalue_3_temp);
		double val_3 = abs(eigenvalue_2_temp - eigenvalue_3_temp);
		if (val_1 <= val_2 && val_1 <= val_3)
		{
			if (abs(eigenvalue_1_temp) > abs(eigenvalue_2_temp))
			{
				eigenvalue_1 = eigenvalue_1_temp;
				eigenvalue_2 = eigenvalue_2_temp;
				eigenvalue_3 = eigenvalue_3_temp;
				int index_1 = 0, index_2 = 1, index_3 = 2;
				eigenvector_1(0) = es.eigenvectors()(0, index_1).real();
				eigenvector_1(1) = es.eigenvectors()(1, index_1).real();
				eigenvector_1(2) = es.eigenvectors()(2, index_1).real();
				eigenvector_2(0) = es.eigenvectors()(0, index_2).real();
				eigenvector_2(1) = es.eigenvectors()(1, index_2).real();
				eigenvector_2(2) = es.eigenvectors()(2, index_2).real();
				eigenvector_3(0) = es.eigenvectors()(0, index_3).real();
				eigenvector_3(1) = es.eigenvectors()(1, index_3).real();
				eigenvector_3(2) = es.eigenvectors()(2, index_3).real();
			}
			else
			{
				eigenvalue_1 = eigenvalue_2_temp;
				eigenvalue_2 = eigenvalue_1_temp;
				eigenvalue_3 = eigenvalue_3_temp;
				int index_1 = 1, index_2 = 0, index_3 = 2;
				eigenvector_1(0) = es.eigenvectors()(0, index_1).real();
				eigenvector_1(1) = es.eigenvectors()(1, index_1).real();
				eigenvector_1(2) = es.eigenvectors()(2, index_1).real();
				eigenvector_2(0) = es.eigenvectors()(0, index_2).real();
				eigenvector_2(1) = es.eigenvectors()(1, index_2).real();
				eigenvector_2(2) = es.eigenvectors()(2, index_2).real();
				eigenvector_3(0) = es.eigenvectors()(0, index_3).real();
				eigenvector_3(1) = es.eigenvectors()(1, index_3).real();
				eigenvector_3(2) = es.eigenvectors()(2, index_3).real();
			}
		}
		else if (val_2 <= val_1 && val_2 <= val_3)
		{
			if (abs(eigenvalue_1_temp) > abs(eigenvalue_3_temp))
			{
				eigenvalue_1 = eigenvalue_1_temp;
				eigenvalue_2 = eigenvalue_3_temp;
				eigenvalue_3 = eigenvalue_2_temp;
				int index_1 = 0, index_2 = 2, index_3 = 1;
				eigenvector_1(0) = es.eigenvectors()(0, index_1).real();
				eigenvector_1(1) = es.eigenvectors()(1, index_1).real();
				eigenvector_1(2) = es.eigenvectors()(2, index_1).real();
				eigenvector_2(0) = es.eigenvectors()(0, index_2).real();
				eigenvector_2(1) = es.eigenvectors()(1, index_2).real();
				eigenvector_2(2) = es.eigenvectors()(2, index_2).real();
				eigenvector_3(0) = es.eigenvectors()(0, index_3).real();
				eigenvector_3(1) = es.eigenvectors()(1, index_3).real();
				eigenvector_3(2) = es.eigenvectors()(2, index_3).real();
			}
			else
			{
				eigenvalue_1 = eigenvalue_3_temp;
				eigenvalue_2 = eigenvalue_1_temp;
				eigenvalue_3 = eigenvalue_2_temp;
				int index_1 = 2, index_2 = 0, index_3 = 1;
				eigenvector_1(0) = es.eigenvectors()(0, index_1).real();
				eigenvector_1(1) = es.eigenvectors()(1, index_1).real();
				eigenvector_1(2) = es.eigenvectors()(2, index_1).real();
				eigenvector_2(0) = es.eigenvectors()(0, index_2).real();
				eigenvector_2(1) = es.eigenvectors()(1, index_2).real();
				eigenvector_2(2) = es.eigenvectors()(2, index_2).real();
				eigenvector_3(0) = es.eigenvectors()(0, index_3).real();
				eigenvector_3(1) = es.eigenvectors()(1, index_3).real();
				eigenvector_3(2) = es.eigenvectors()(2, index_3).real();
			}
		}
		else
		{
			if (abs(eigenvalue_2_temp) > abs(eigenvalue_3_temp))
			{
				eigenvalue_1 = eigenvalue_2_temp;
				eigenvalue_2 = eigenvalue_3_temp;
				eigenvalue_3 = eigenvalue_1_temp;
				int index_1 = 1, index_2 = 2, index_3 = 0;
				eigenvector_1(0) = es.eigenvectors()(0, index_1).real();
				eigenvector_1(1) = es.eigenvectors()(1, index_1).real();
				eigenvector_1(2) = es.eigenvectors()(2, index_1).real();
				eigenvector_2(0) = es.eigenvectors()(0, index_2).real();
				eigenvector_2(1) = es.eigenvectors()(1, index_2).real();
				eigenvector_2(2) = es.eigenvectors()(2, index_2).real();
				eigenvector_3(0) = es.eigenvectors()(0, index_3).real();
				eigenvector_3(1) = es.eigenvectors()(1, index_3).real();
				eigenvector_3(2) = es.eigenvectors()(2, index_3).real();
			}
			else
			{
				eigenvalue_1 = eigenvalue_3_temp;
				eigenvalue_2 = eigenvalue_2_temp;
				eigenvalue_3 = eigenvalue_1_temp;
				int index_1 = 2, index_2 = 1, index_3 = 0;
				eigenvector_1(0) = es.eigenvectors()(0, index_1).real();
				eigenvector_1(1) = es.eigenvectors()(1, index_1).real();
				eigenvector_1(2) = es.eigenvectors()(2, index_1).real();
				eigenvector_2(0) = es.eigenvectors()(0, index_2).real();
				eigenvector_2(1) = es.eigenvectors()(1, index_2).real();
				eigenvector_2(2) = es.eigenvectors()(2, index_2).real();
				eigenvector_3(0) = es.eigenvectors()(0, index_3).real();
				eigenvector_3(1) = es.eigenvectors()(1, index_3).real();
				eigenvector_3(2) = es.eigenvectors()(2, index_3).real();
			}
		}
	}
	if (eigenvalue_1_temp <= 0 && eigenvalue_2_temp <= 0 && eigenvalue_3_temp <= 0)
	{
	}
	if (eigenvector_3(2) > 0)
		eigenvector_1 = eigenvector_2.cross(eigenvector_3);
	else
	{
		eigenvector_3 = -1.0 * eigenvector_3;
		eigenvector_1 = eigenvector_2.cross(eigenvector_3);
	}


	double Ell_flag_a = sqrt((abs(eigenvalue_1) - abs(eigenvalue_2)) / (abs(eigenvalue_1) + abs(eigenvalue_3)));
	double Ell_flag_b = sqrt((abs(eigenvalue_2) + abs(eigenvalue_3)) / (abs(eigenvalue_1) + abs(eigenvalue_3)));

	Eigen::Vector3d Ell_N_1;
	Eigen::Vector3d Ell_N_2;
	Ell_N_1(0) = eigenvector_1(0) * Ell_flag_a - eigenvector_3(0) * Ell_flag_b;
	Ell_N_1(1) = eigenvector_1(1) * Ell_flag_a - eigenvector_3(1) * Ell_flag_b;
	Ell_N_1(2) = eigenvector_1(2) * Ell_flag_a - eigenvector_3(2) * Ell_flag_b;
	Ell_N_2(0) = -eigenvector_1(0) * Ell_flag_a - eigenvector_3(0) * Ell_flag_b;
	Ell_N_2(1) = -eigenvector_1(1) * Ell_flag_a - eigenvector_3(1) * Ell_flag_b;
	Ell_N_2(2) = -eigenvector_1(2) * Ell_flag_a - eigenvector_3(2) * Ell_flag_b;
	double M_1 = sqrt(Ell_N_1(0) * Ell_N_1(0) + Ell_N_1(1) * Ell_N_1(1) + Ell_N_1(2) * Ell_N_1(2));
	double M_2 = sqrt(Ell_N_2(0) * Ell_N_2(0) + Ell_N_2(1) * Ell_N_2(1) + Ell_N_2(2) * Ell_N_2(2));
	Ell_N_1 = Ell_N_1 / M_1;
	Ell_N_2 = Ell_N_2 / M_2;
	ellipse_view_pose as;
	as.view_N_1[0] = Ell_N_1(0);
	as.view_N_1[1] = Ell_N_1(1);
	as.view_N_1[2] = Ell_N_1(2);
	as.view_N_2[0] = Ell_N_2(0);
	as.view_N_2[1] = Ell_N_2(1);
	as.view_N_2[2] = Ell_N_2(2);
	return as;
}
void SFM_module::update_points_in_point(std::vector<std::vector<Coded_detect_inf>> keypoints_in_queue,
	std::vector<Coded_detect_inf> keypoints_now, cv::Mat K, cv::Mat Dis,
	std::vector<cv::Mat> R_in_queue, std::vector<cv::Mat> T_in_queue, cv::Mat R_now, cv::Mat T_now,
	std::vector<sfm_3d_group>& P_world, double& update_number, int& new_number)
{
	update_number = 1;
	new_number = 0;
	//遍历和之前图像的交集，通过F矩阵检查才认为质量通过，可以计算。
	for (int ii = 0; ii < keypoints_in_queue.size(); ii++)
	{
		std::vector<cv::Point2d> P1_ori, P2_ori;
		std::vector<int> Code_ori;
		cv::Mat R_ori, T_ori;	
		cv::Mat mask_ori;
		get_matched_points(keypoints_in_queue[ii],
			keypoints_now, P1_ori, P2_ori, Code_ori);
		if (P1_ori.size() <= 5)
		{
			continue;
		}
		double focal_length = 0.5 * (K.at<double>(0) + K.at<double>(4));
		cv::Point2d principle_point(K.at<double>(2), K.at<double>(5));

		cv::Mat E = cv::findEssentialMat(P1_ori, P2_ori, focal_length, principle_point, cv::RANSAC, 0.999, 3.0, mask_ori);
		if (E.empty())
		{
			continue;
		}
		double feasible_count = countNonZero(mask_ori);
		if (feasible_count <= 5 || (feasible_count / P1_ori.size()) < 0.6)
		{
			continue;
		}

		int pass_count = recoverPose(E, P1_ori, P2_ori, R_ori, T_ori, focal_length, principle_point, mask_ori);
		if (((double)pass_count) / feasible_count < 0.6)
		{
			continue;
		}
		update_number = update_number > ((double)pass_count) / feasible_count ? ((double)pass_count) / feasible_count : update_number;
		maskout_points(P1_ori, P2_ori, Code_ori, mask_ori);
		std::vector<sfm_3d_group> xyz_ori = reconstruct_3d_continue(K, Dis, R_in_queue[ii], T_in_queue[ii], R_now, T_now,
			P1_ori, P2_ori, Code_ori);
		for (int jj = 0; jj < xyz_ori.size(); jj++)
		{
			bool exist_code=false;
			for (int pp = 0; pp < P_world.size(); pp++)
			{
				if (xyz_ori[jj].code == P_world[pp].code)
				{
					exist_code = true;
					break;
				}
			}
			if (!exist_code)
			{
				new_number++;
				P_world.push_back(xyz_ori[jj]);
			}
		}
	}
}
void SFM_module::remove_calculated_index_from_pair(std::vector<std::vector<Coded_detect_inf>>& P_serial, std::vector<std::vector<std::vector<cv::Point>>>& PC_serial, std::vector<int>& P_serial_index,
	match_pair_group matches)
{
	std::vector<std::vector<Coded_detect_inf>> P_serial_copy = P_serial;
	std::vector<std::vector<std::vector<cv::Point>>> PC_serial_copy = PC_serial;
	std::vector<int> P_serial_index_copy = P_serial_index;
	P_serial.clear();
	PC_serial.clear();
	P_serial_index.clear();
	for (int ii = 0; ii < P_serial_copy.size(); ii++)
	{
		if (ii != matches.position_x && ii != matches.position_y)
		{
			P_serial.push_back(P_serial_copy[ii]);
			PC_serial.push_back(PC_serial_copy[ii]);
			P_serial_index.push_back(P_serial_index_copy[ii]);
		}
	}
}
void SFM_module::remove_calculated_index_from_index(std::vector<std::vector<Coded_detect_inf>>& P_serial, std::vector<std::vector<std::vector<cv::Point>>>& PC_serial, std::vector<int>& P_serial_index,
	int match)
{
	std::vector<std::vector<Coded_detect_inf>> P_serial_copy = P_serial;
	std::vector<std::vector<std::vector<cv::Point>>> PC_serial_copy = PC_serial;
	std::vector<int> P_serial_index_copy = P_serial_index;
	P_serial.clear();
	PC_serial.clear();
	P_serial_index.clear();
	for (int ii = 0; ii < P_serial_copy.size(); ii++)
	{
		if (ii != match)
		{
			P_serial.push_back(P_serial_copy[ii]);
			PC_serial.push_back(PC_serial_copy[ii]);
			P_serial_index.push_back(P_serial_index_copy[ii]);
		}
	}
}
std::vector<sfm_3d_group> SFM_module::reconstruct_3d_first(cv::Mat K, cv::Mat R, cv::Mat T, std::vector<cv::Point2d> P1, std::vector<cv::Point2d> P2, std::vector<int> Code)
{
	//两个相机的投影矩阵[R T]，triangulatePoints只支持float型
	cv::Mat proj1(3, 4, CV_64FC1);
	cv::Mat proj2(3, 4, CV_64FC1);

	proj1(cv::Range(0, 3), cv::Range(0, 3)) = cv::Mat::eye(3, 3, CV_64FC1);
	proj1.col(3) = cv::Mat::zeros(3, 1, CV_64FC1);

	R.convertTo(proj2(cv::Range(0, 3), cv::Range(0, 3)), CV_64FC1);
	T.convertTo(proj2.col(3), CV_64FC1);

	proj1 = K * proj1;
	proj2 = K * proj2;

	std::vector<sfm_3d_group> d3_group;
	cv::Mat structure;
	triangulatePoints(proj1, proj2, P1, P2, structure);
	for (int ii = 0; ii < P1.size(); ii++)
	{
		double x = structure.at<double>(0, ii) / structure.at<double>(3, ii);
		double y = structure.at<double>(1, ii) / structure.at<double>(3, ii);
		double z = structure.at<double>(2, ii) / structure.at<double>(3, ii);
		d3_group.push_back(sfm_3d_group(Code[ii], 1, x, y, z, 0));
	}
	return d3_group;
}
std::vector<sfm_3d_group> SFM_module::reconstruct_3d_continue(cv::Mat K, cv::Mat Dis, cv::Mat R1, cv::Mat T1, cv::Mat R2, cv::Mat T2,
	std::vector<cv::Point2d> P1, std::vector<cv::Point2d> P2, std::vector<int> Code)
{
	cv::Mat proj1(3, 4, CV_64FC1);
	cv::Mat proj2(3, 4, CV_64FC1);


	R1.convertTo(proj1(cv::Range(0, 3), cv::Range(0, 3)), CV_64FC1);
	T1.convertTo(proj1.col(3), CV_64FC1);
	R2.convertTo(proj2(cv::Range(0, 3), cv::Range(0, 3)), CV_64FC1);
	T2.convertTo(proj2.col(3), CV_64FC1);

	proj1 = K * proj1;
	proj2 = K * proj2;
	cv::undistortPoints(P1, P1, K, Dis, K);
	cv::undistortPoints(P2, P2, K, Dis, K);
	std::vector<sfm_3d_group> d3_group;
	cv::Mat structure;
	triangulatePoints(proj1, proj2, P1, P2, structure);
	for (int ii = 0; ii < P1.size(); ii++)
	{
		double x = structure.at<double>(0, ii) / structure.at<double>(3, ii);
		double y = structure.at<double>(1, ii) / structure.at<double>(3, ii);
		double z = structure.at<double>(2, ii) / structure.at<double>(3, ii);
		d3_group.push_back(sfm_3d_group(Code[ii], 1, x, y, z, 0));
	}
	return d3_group;
}
void SFM_module::maskout_points(std::vector<cv::Point2d> &P_1, std::vector<cv::Point2d> &P_2, std::vector<int>& Code, cv::Mat mask)
{
	std::vector<cv::Point2d> P_1_copy = P_1, P_2_copy = P_2;
	std::vector<int> Code_copy = Code;
	P_1.clear();
	P_2.clear();
	Code.clear();

	for (int ii = 0; ii < mask.rows; ++ii)
	{
		if (mask.at<uchar>(ii) > 0)
		{
			P_1.push_back(P_1_copy[ii]);
			P_2.push_back(P_2_copy[ii]);
			Code.push_back(Code_copy[ii]);
		}
	}
}

void SFM_module::get_matched_points(std::vector<Coded_detect_inf> P1, std::vector<Coded_detect_inf>P2,
	std::vector<cv::Point2d>& out_P1,std::vector<cv::Point2d>& out_P2, std::vector<int>& Code)
{
	out_P1.clear();
	out_P2.clear();
	Code.clear();
	for (int ii = 0; ii < P1.size(); ii++)
	{
		for (int jj = 0; jj < P2.size(); jj++)
		{
			if (P1[ii].code_num == P2[jj].code_num)
			{
				out_P1.push_back(cv::Point2d(P1[ii].x, P1[ii].y));
				out_P2.push_back(cv::Point2d(P2[jj].x, P2[jj].y));
				Code.push_back(P1[ii].code_num);
			}
		}
	}
}

int SFM_module::cal_3d_transform(cv::Mat& K, cv::Mat& Dis, std::vector<sfm_3d_group> D3_point, std::vector<Coded_detect_inf> D2_point,
	std::vector<cv::Point3d> &out_3d, std::vector<cv::Point2d> &out_2d, std::vector<int> &code, cv::Mat& R, cv::Mat& T, cv::Mat& mask)
{
	out_3d.clear();
	out_2d.clear();
	code.clear();
	cv::Mat R_temp;
	for (int ii = 0; ii < D3_point.size(); ii++)
	{
		for (int jj = 0; jj < D2_point.size(); jj++)
		{
			if (D3_point[ii].code == D2_point[jj].code_num)
			{
				out_3d.push_back(cv::Point3d(D3_point[ii].x, D3_point[ii].y, D3_point[ii].z));
				out_2d.push_back(cv::Point2d(D2_point[jj].x, D2_point[jj].y));
				code.push_back(D3_point[ii].code);
			}
		}
	}
	//cv::solvePnP(out_3d, out_2d, K, Dis, R_temp, T, false);
	//Rodrigues(R_temp, R);
	//return out_2d.size();

	cv::solvePnPRansac(out_3d, out_2d, K, Dis, R_temp, T, false, 1000, SFM_setting_window->ui.pnp_error_doubleSpinBox->value(), 0.99, mask, cv::SOLVEPNP_ITERATIVE);
	Rodrigues(R_temp, R);
	double feasible_count = mask.rows;
	if (feasible_count > 4)
	{
		Eigen::Matrix<double, 3, 3> matrix_temp;
		cv::cv2eigen(R, matrix_temp);
		Eigen::Quaterniond R_temp_s = Eigen::Quaterniond(matrix_temp);
		R_temp_s.normalize();
		Eigen::Matrix<double, 3, 1> matrix_temp2;
		cv::cv2eigen(T, matrix_temp2);
		Eigen::Vector3d T_temp_s = Eigen::Vector3d(matrix_temp2);
		double camera_k[5];
		double camera_dis_k[6];
		double camera_dis_p[2];
		double camera_dis_t[4];
		camera_k[0] = K.at<double>(0, 0); camera_k[1] = K.at<double>(1, 1); camera_k[2] = K.at<double>(0, 2); camera_k[3] = K.at<double>(1, 2); camera_k[4] = K.at<double>(0, 1);
		camera_dis_k[0] = Dis.at<double>(0, 0);
		camera_dis_k[1] = Dis.at<double>(0, 1);
		camera_dis_k[2] = Dis.at<double>(0, 4);
		camera_dis_k[3] = Dis.at<double>(0, 5);
		camera_dis_k[4] = Dis.at<double>(0, 6);
		camera_dis_k[5] = Dis.at<double>(0, 7);
		camera_dis_p[0] = Dis.at<double>(0, 2);
		camera_dis_p[1] = Dis.at<double>(0, 3);
		camera_dis_t[0] = Dis.at<double>(0, 8);
		camera_dis_t[2] = Dis.at<double>(0, 9);
		camera_dis_t[1] = Dis.at<double>(0, 10);
		camera_dis_t[3] = Dis.at<double>(0, 11);
		SFM_optimize::calculate_RT(out_3d, out_2d, R_temp_s, T_temp_s, camera_k, camera_dis_k, camera_dis_p, camera_dis_t
			, SFM_setting_window->ui.max_iter_spinBox->value()
			, ((double)SFM_setting_window->ui.stop_spinBox->value()) * pow(10, SFM_setting_window->ui.stop_spinBox_2->value())
			, SFM_setting_window->ui.thread_spinBox->value()
			, SFM_setting_window->ui.timeout_spinBox->value()
			, SFM_setting_window->ui.loss_comboBox->currentIndex(), SFM_setting_window->ui.pnp_error_doubleSpinBox->value());
		cv::eigen2cv(R_temp_s.matrix(), R);
		cv::eigen2cv(T_temp_s, T);
	}
	return feasible_count;
}
int SFM_module::cal_transform(cv::Mat& K, std::vector<cv::Point2d>& p1, std::vector<cv::Point2d>& p2, cv::Mat& R, cv::Mat& T, cv::Mat& mask, int mode)
{
	bool has_found_from_E = false;
	bool has_found_from_H = false;
	cv::Mat R_from_E; cv::Mat T_from_E;
	cv::Mat R_from_H; cv::Mat T_from_H;
	cv::Mat mask_from_E; cv::Mat mask_from_H;
	int pass_count = 0;
	int pass_count_from_E = 0;
	int pass_count_from_H = 0;
	//根据内参矩阵获取相机的焦距和光心坐标（主点坐标）
	double focal_length = 0.5 * (K.at<double>(0) + K.at<double>(4));
	cv::Point2d principle_point(K.at<double>(2), K.at<double>(5));

	if (mode == 0 || mode == 2)
	{
		//根据匹配点求取本征矩阵，使用RANSAC，进一步排除失配点
		cv::Mat E = findEssentialMat(p1, p2, focal_length, principle_point, cv::RANSAC, 0.999, 3, mask_from_E);
		if (!E.empty())
		{
			double feasible_count = countNonZero(mask_from_E);
			//对于RANSAC而言，outlier数量大于50%时，结果是不可靠的
			if (!(feasible_count <= 5 || (feasible_count / p1.size()) < 0.6))
			{
				//分解本征矩阵，获取相对变换
				pass_count_from_E = recoverPose(E, p1, p2, R_from_E, T_from_E, focal_length, principle_point, mask_from_E);
				//同时位于两个相机前方的点的数量要足够大
				if (!(((double)pass_count_from_E) / feasible_count < 0.7))
				{
					has_found_from_E = true;
				}
			}
		}
	}
	if (mode == 1 || mode == 2)
	{
		cv::Mat H = cv::findHomography(p1, p2, cv::RANSAC, 3, mask_from_H, 2000, 0.995);
		if (!H.empty())
		{
			pass_count_from_H = countNonZero(mask_from_H);
			if (!(pass_count_from_H <= 5 || (pass_count_from_H / p1.size()) < 0.6))
			{
				std::vector<cv::Mat> Rs_decomp, ts_decomp, normals_decomp;
				std::vector<cv::Point2f> befor_p;
				std::vector<cv::Point2f> after_p;
				for (int ii = 0; ii < p1.size(); ++ii)
				{
					befor_p.push_back(cv::Point2f(p1[ii].x, p1[ii].y));
					after_p.push_back(cv::Point2f(p2[ii].x, p2[ii].y));
				}
				int solutions = cv::decomposeHomographyMat(H, K, Rs_decomp, ts_decomp, normals_decomp);
				std::vector<int>  sol;
				filterHomographyDecompByVisibleRefpoints(Rs_decomp, normals_decomp, befor_p, after_p, sol, mask_from_H);
				//for (int ii = 0; ii < Rs_decomp.size(); ii++)
				//{
				//	cv::Mat proj1_L(3, 4, CV_64FC1);
				//	cv::Mat proj2_L(3, 4, CV_64FC1);
				//	cv::Mat proj1_R(3, 4, CV_64FC1);
				//	cv::Mat proj2_R(3, 4, CV_64FC1);

				//	proj1_L(cv::Range(0, 3), cv::Range(0, 3)) = cv::Mat::eye(3, 3, CV_64FC1);
				//	proj1_L.col(3) = cv::Mat::zeros(3, 1, CV_64FC1);
				//	proj1_R(cv::Range(0, 3), cv::Range(0, 3)) = cv::Mat::eye(3, 3, CV_64FC1);
				//	proj1_R.col(3) = cv::Mat::zeros(3, 1, CV_64FC1);

				//	cv::Mat mat_R_L = Rs_decomp[ii];
				//	cv::Mat mat_T_L = ts_decomp[ii];
				//	mat_R_L.convertTo(proj2_L(cv::Range(0, 3), cv::Range(0, 3)), CV_64FC1);
				//	mat_T_L.convertTo(proj2_L.col(3), CV_64FC1);

				//	cv::Mat mat_R_R = Rs_decomp[ii].inv();
				//	cv::Mat mat_T_R = Rs_decomp[ii].inv() * ts_decomp[ii];
				//	mat_R_R.convertTo(proj2_R(cv::Range(0, 3), cv::Range(0, 3)), CV_64FC1);
				//	mat_T_R.convertTo(proj2_R.col(3), CV_64FC1);

				//	proj1_L = K * proj1_L;
				//	proj2_L = K * proj2_L;
				//	proj1_R = K * proj1_R;
				//	proj2_R = K * proj2_R;

				//	cv::Mat structure_L;
				//	std::vector<cv::Point3d> xyz_now_L;
				//	cv::Mat structure_R;
				//	std::vector<cv::Point3d> xyz_now_R;
				//	triangulatePoints(proj1_L, proj2_L, p1, p2, structure_L);
				//	triangulatePoints(proj1_R, proj2_R, p2, p1, structure_R);
				//	for (int jj = 0; jj < p1.size(); jj++)
				//	{
				//		xyz_now_L.push_back(cv::Point3d(structure_L.at<double>(0, jj) / structure_L.at<double>(3, jj), structure_L.at<double>(1, jj) / structure_L.at<double>(3, jj), structure_L.at<double>(2, jj) / structure_L.at<double>(3, jj)));
				//		xyz_now_R.push_back(cv::Point3d(structure_R.at<double>(0, jj) / structure_R.at<double>(3, jj), structure_R.at<double>(1, jj) / structure_R.at<double>(3, jj), structure_R.at<double>(2, jj) / structure_R.at<double>(3, jj)));
				//	}
				//}
				if (sol.size() != 0)
				{
					R_from_H = Rs_decomp[sol[0]];
					T_from_H = ts_decomp[sol[0]];
					mask_from_E = mask_from_H;
					has_found_from_H = true;
				}
			}
		}
	}
	if ((has_found_from_E && has_found_from_H && pass_count_from_E > pass_count_from_H)|| (has_found_from_E && !has_found_from_H))
	{
		pass_count = pass_count_from_E;
		R = R_from_E;
		T = T_from_E;
		mask = mask_from_E;
	}
	else if(has_found_from_H)
	{
		pass_count = pass_count_from_H;
		R = R_from_H;
		T = T_from_H;
		mask = mask_from_H;
	}
	return pass_count;
}

int SFM_module::find_same_coded_point_to_3d(std::vector<sfm_3d_group> P_world, std::vector<Coded_detect_inf> P)
{
	int same_number = 0;
	for (int ii = 0; ii < P_world.size(); ii++)
	{
		for (int jj = 0; jj < P.size(); jj++)
		{
			if (P_world[ii].code == P[jj].code_num)
			{
				same_number++;
			}
		}
	}
	return same_number;
}
int SFM_module::find_same_coded_point(std::vector<Coded_detect_inf> P1, std::vector<Coded_detect_inf>P2)
{
	int same_number = 0;
	for (int ii = 0; ii < P1.size(); ii++)
	{
		for (int jj = 0; jj < P2.size(); jj++)
		{
			if (P1[ii].code_num == P2[jj].code_num)
			{
				same_number++;
			}
		}
	}
	return same_number;
}

void SFM_module::read_KeyPoint_file(QString path_read)
{
	in_reading_keypoint = true;
	FILE* fp;
	size_t pointer_in_file = 0;
	fp = fopen(path_read.toLocal8Bit().data(), "rb");
	if (NULL == fp)
	{
		QMessageBox::warning(this, tr("Read KeyPoint File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
		clear_all_data();
		in_reading_keypoint = false;
		update_show_view();
		update_table_view();
		ui.open_path_lineEdit->setText(tr(""));
		fclose(fp);
		return;
	}
	unsigned char read_file_mode = -1;
	size_t return_size = fread(&read_file_mode, sizeof(unsigned char), 1, fp);
	if (return_size != 1)
	{
		QMessageBox::warning(this, tr("Read KeyPoint File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
		clear_all_data();
		update_show_view();
		in_reading_keypoint = false;
		update_table_view();
		ui.open_path_lineEdit->setText(tr(""));
		fclose(fp);
		return;
	}
	if (read_file_mode != 31)
	{
		QMessageBox::warning(this, tr("Read KeyPoint File"), tr("Storage content not SFM-KeyPoint data!"), QMessageBox::Ok, QMessageBox::Ok);
		clear_all_data();
		update_show_view();
		in_reading_keypoint = false;
		update_table_view();
		ui.open_path_lineEdit->setText(tr(""));
		fclose(fp);
		return;
	}
	if (feof(fp))
	{
		QMessageBox::warning(this, tr("Read KeyPoint File"), tr("Read Successfully, but no data"), QMessageBox::Ok, QMessageBox::Ok);
		clear_all_data();
		update_show_view();
		in_reading_keypoint = false;
		update_table_view();
		ui.open_path_lineEdit->setText(tr(""));
		fclose(fp);
		return;
	}
	while (!feof(fp))
	{
		int image_index;
		return_size = fread(&image_index, sizeof(int), 1, fp);
		if (return_size != 1)
		{
			if (feof(fp))
			{
				break;
			}
			QMessageBox::warning(this, tr("Read KeyPoint File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			clear_all_data();
			update_show_view();
			in_reading_keypoint = false;
			update_table_view();
			ui.open_path_lineEdit->setText(tr(""));
			fclose(fp);
			return;
		}
		long long string_path_length;
		return_size = fread(&string_path_length, sizeof(long long), 1, fp);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Read KeyPoint File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			clear_all_data();
			update_show_view();
			in_reading_keypoint = false;
			update_table_view();
			ui.open_path_lineEdit->setText(tr(""));
			fclose(fp);
			return;
		}
		std::string string_path;
		string_path.resize(string_path_length);
		return_size = fread(&string_path[0], sizeof(char), string_path_length, fp);
		if (return_size != string_path_length)
		{
			QMessageBox::warning(this, tr("Read KeyPoint File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			clear_all_data();
			update_show_view();
			in_reading_keypoint = false;
			update_table_view();
			ui.open_path_lineEdit->setText(tr(""));
			fclose(fp);
			return;
		}
		Image_serial_name.push_back(QString(QString::fromLocal8Bit(string_path.data())));
		KeyPoint_Enable_serial.push_back(true);
		QDir file_temp = QDir(QString(QString::fromLocal8Bit(string_path.data())));
		file_temp.cdUp();
		ui.open_path_lineEdit->setText(file_temp.absolutePath());
		int image_width = -1, image_height = -1;
		return_size = fread(&image_width, sizeof(int), 1, fp);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Read KeyPoint File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			clear_all_data();
			update_show_view();
			in_reading_keypoint = false;
			update_table_view();
			ui.open_path_lineEdit->setText(tr(""));
			fclose(fp);
			return;
		}
		return_size = fread(&image_height, sizeof(int), 1, fp);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Read KeyPoint File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			clear_all_data();
			update_show_view();
			in_reading_keypoint = false;
			update_table_view();
			ui.open_path_lineEdit->setText(tr(""));
			fclose(fp);
			return;
		}
		ImageWidth_serial.push_back(image_width);
		ImageHeight_serial.push_back(image_height);
		int corn_size = -1;
		return_size = fread(&corn_size, sizeof(int), 1, fp);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Read KeyPoint File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			clear_all_data();
			update_show_view();
			in_reading_keypoint = false;
			update_table_view();
			ui.open_path_lineEdit->setText(tr(""));
			fclose(fp);
			return;
		}
		std::vector<Coded_detect_inf> kp_temp;
		std::vector<std::vector<cv::Point>> contours_temp;
		for (int pp = 0; pp < corn_size; pp++)
		{
			double temp_x, temp_y, temp_ra, temp_rb, temp_theta, temp_err;
			double temp_N1, temp_N2, temp_N3, temp_N4, temp_N5, temp_N6;
			int temp_code,temp_contours_num;
			return_size = fread(&temp_x, sizeof(double), 1, fp);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Read KeyPoint File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				clear_all_data();
				update_show_view();
				in_reading_keypoint = false;
				update_table_view();
				ui.open_path_lineEdit->setText(tr(""));
				fclose(fp);
				return;
			}
			return_size = fread(&temp_y, sizeof(double), 1, fp);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Read KeyPoint File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				clear_all_data();
				update_show_view();
				in_reading_keypoint = false;
				update_table_view();
				ui.open_path_lineEdit->setText(tr(""));
				fclose(fp);
				return;
			}
			return_size = fread(&temp_ra, sizeof(double), 1, fp);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Read KeyPoint File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				clear_all_data();
				update_show_view();
				in_reading_keypoint = false;
				update_table_view();
				ui.open_path_lineEdit->setText(tr(""));
				fclose(fp);
				return;
			}
			return_size = fread(&temp_rb, sizeof(double), 1, fp);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Read KeyPoint File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				clear_all_data();
				update_show_view();
				in_reading_keypoint = false;
				update_table_view();
				ui.open_path_lineEdit->setText(tr(""));
				fclose(fp);
				return;
			}
			return_size = fread(&temp_theta, sizeof(double), 1, fp);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Read KeyPoint File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				clear_all_data();
				update_show_view();
				in_reading_keypoint = false;
				update_table_view();
				ui.open_path_lineEdit->setText(tr(""));
				fclose(fp);
				return;
			}
			return_size = fread(&temp_code, sizeof(int), 1, fp);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Read KeyPoint File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				clear_all_data();
				update_show_view();
				in_reading_keypoint = false;
				update_table_view();
				ui.open_path_lineEdit->setText(tr(""));
				fclose(fp);
				return;
			}
			return_size = fread(&temp_err, sizeof(double), 1, fp);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Read KeyPoint File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				clear_all_data();
				update_show_view();
				in_reading_keypoint = false;
				update_table_view();
				ui.open_path_lineEdit->setText(tr(""));
				fclose(fp);
				return;
			}
			return_size = fread(&temp_N1, sizeof(double), 1, fp);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Read KeyPoint File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				clear_all_data();
				update_show_view();
				in_reading_keypoint = false;
				update_table_view();
				ui.open_path_lineEdit->setText(tr(""));
				fclose(fp);
				return;
			}
			return_size = fread(&temp_N2, sizeof(double), 1, fp);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Read KeyPoint File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				clear_all_data();
				update_show_view();
				in_reading_keypoint = false;
				update_table_view();
				ui.open_path_lineEdit->setText(tr(""));
				fclose(fp);
				return;
			}
			return_size = fread(&temp_N3, sizeof(double), 1, fp);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Read KeyPoint File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				clear_all_data();
				update_show_view();
				in_reading_keypoint = false;
				update_table_view();
				ui.open_path_lineEdit->setText(tr(""));
				fclose(fp);
				return;
			}
			return_size = fread(&temp_N4, sizeof(double), 1, fp);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Read KeyPoint File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				clear_all_data();
				update_show_view();
				in_reading_keypoint = false;
				update_table_view();
				ui.open_path_lineEdit->setText(tr(""));
				fclose(fp);
				return;
			}
			return_size = fread(&temp_N5, sizeof(double), 1, fp);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Read KeyPoint File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				clear_all_data();
				update_show_view();
				in_reading_keypoint = false;
				update_table_view();
				ui.open_path_lineEdit->setText(tr(""));
				fclose(fp);
				return;
			}
			return_size = fread(&temp_N6, sizeof(double), 1, fp);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Read KeyPoint File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				clear_all_data();
				update_show_view();
				in_reading_keypoint = false;
				update_table_view();
				ui.open_path_lineEdit->setText(tr(""));
				fclose(fp);
				return;
			}
			ellipse_view_pose temp_pose;
			temp_pose.view_N_1[0] = temp_N1;
			temp_pose.view_N_1[1] = temp_N2;
			temp_pose.view_N_1[2] = temp_N3;
			temp_pose.view_N_2[0] = temp_N4;
			temp_pose.view_N_2[1] = temp_N5;
			temp_pose.view_N_2[2] = temp_N6;
			Coded_detect_inf points_now_temp;
			points_now_temp.x = temp_x;
			points_now_temp.y = temp_y;
			points_now_temp.r_a = temp_ra;
			points_now_temp.r_b = temp_rb;
			points_now_temp.angle_in_pi = temp_theta;
			points_now_temp.code_num = temp_code;
			points_now_temp.fit_err = temp_err;
			points_now_temp.pose = temp_pose;
			kp_temp.push_back(points_now_temp);
			std::vector<cv::Point> contours_now;
			return_size = fread(&temp_contours_num, sizeof(int), 1, fp);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Read KeyPoint File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				clear_all_data();
				update_show_view();
				in_reading_keypoint = false;
				update_table_view();
				ui.open_path_lineEdit->setText(tr(""));
				fclose(fp);
				return;
			}
			for (int qq = 0; qq < temp_contours_num; qq++)
			{
				int temp_x_contour, temp_y_contour;
				return_size = fread(&temp_x_contour, sizeof(int), 1, fp);
				if (return_size != 1)
				{
					QMessageBox::warning(this, tr("Read KeyPoint File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
					clear_all_data();
					update_show_view();
					in_reading_keypoint = false;
					update_table_view();
					ui.open_path_lineEdit->setText(tr(""));
					fclose(fp);
					return;
				}
				return_size = fread(&temp_y_contour, sizeof(int), 1, fp);
				if (return_size != 1)
				{
					QMessageBox::warning(this, tr("Read KeyPoint File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
					clear_all_data();
					update_show_view();
					in_reading_keypoint = false;
					update_table_view();
					ui.open_path_lineEdit->setText(tr(""));
					fclose(fp);
					return;
				}
				contours_now.push_back(cv::Point(temp_x_contour, temp_y_contour));
			}
			contours_temp.push_back(contours_now);
		}
		Result_contours_serial.push_back(contours_temp);
		KeyPoint_code_serial.push_back(kp_temp);
	}

	if (Image_serial_name.size() == 0)
	{
		ui.show_detect_spinBox->setEnabled(false);
		ui.show_detect_horizontalSlider->setEnabled(false);
		ui.show_detect_spinBox->setMinimum(1);
		ui.show_detect_spinBox->setMaximum(1);
		ui.show_detect_horizontalSlider->setMinimum(1);
		ui.show_detect_horizontalSlider->setMaximum(1);
		ui.show_detect_spinBox->setValue(1);
		ui.show_detect_horizontalSlider->setValue(1);
		ui.Detect_view->clear_image();
	}
	else
	{
		ui.show_detect_spinBox->setEnabled(true);
		ui.show_detect_horizontalSlider->setEnabled(true);
		ui.show_detect_spinBox->setMinimum(1);
		ui.show_detect_spinBox->setMaximum(Image_serial_name.size());
		ui.show_detect_horizontalSlider->setMinimum(1);
		ui.show_detect_horizontalSlider->setMaximum(Image_serial_name.size());
		ui.show_detect_spinBox->setValue(1);
		ui.show_detect_horizontalSlider->setValue(1);
	}

	in_reading_keypoint = false;
	need_recal_keypoints = false;
	update_show_view();
	update_table_view();
	QMessageBox::information(this, tr("Read KeyPoint File"), tr("Read Successfully."), QMessageBox::Ok);
	fclose(fp);
}
void SFM_module::read_SFMResult_file(QString path_read)
{
	in_reading_keypoint = true;
	FILE* fp;
	size_t pointer_in_file = 0;
	fp = fopen(path_read.toLocal8Bit().data(), "rb");
	if (NULL == fp)
	{
		QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
		clear_all_data();
		in_reading_keypoint = false;
		update_show_view();
		update_points_view();
		update_table_view();
		ui.open_path_lineEdit->setText(tr(""));
		fclose(fp);
		return;
	}
	unsigned char read_file_mode = -1;
	size_t return_size = fread(&read_file_mode, sizeof(unsigned char), 1, fp);
	if (return_size != 1)
	{
		QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
		clear_all_data();
		update_show_view();
		in_reading_keypoint = false;
		update_points_view();
		update_table_view();
		ui.open_path_lineEdit->setText(tr(""));
		fclose(fp);
		return;
	}
	if (read_file_mode != 32)
	{
		QMessageBox::warning(this, tr("Read SFM Result File"), tr("Storage content not SFM-KeyPoint data!"), QMessageBox::Ok, QMessageBox::Ok);
		clear_all_data();
		update_show_view();
		in_reading_keypoint = false;
		update_points_view();
		update_table_view();
		ui.open_path_lineEdit->setText(tr(""));
		fclose(fp);
		return;
	}
	double* camera_K = new double[5];
	return_size = fread(camera_K, sizeof(double), 5, fp);
	if (return_size != 5)
	{
		QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Successfully, but no data"), QMessageBox::Ok, QMessageBox::Ok);
		clear_all_data();
		update_show_view();
		in_reading_keypoint = false;
		update_points_view();
		update_table_view();
		ui.open_path_lineEdit->setText(tr(""));
		fclose(fp);
		return;
	}
	Result_Camera_K[0] = camera_K[0];
	Result_Camera_K[1] = camera_K[1];
	Result_Camera_K[2] = camera_K[2];
	Result_Camera_K[3] = camera_K[3];
	Result_Camera_K[4] = camera_K[4];
	delete[] camera_K;

	double* camera_Dis_K = new double[6];
	return_size = fread(camera_Dis_K, sizeof(double), 6, fp);
	if (return_size != 6)
	{
		QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Successfully, but no data"), QMessageBox::Ok, QMessageBox::Ok);
		clear_all_data();
		update_show_view();
		in_reading_keypoint = false;
		update_points_view();
		update_table_view();
		ui.open_path_lineEdit->setText(tr(""));
		fclose(fp);
		return;
	}

	Result_Camera_Dis_K[0] = camera_Dis_K[0];
	Result_Camera_Dis_K[1] = camera_Dis_K[1];
	Result_Camera_Dis_K[2] = camera_Dis_K[2];
	Result_Camera_Dis_K[3] = camera_Dis_K[3];
	Result_Camera_Dis_K[4] = camera_Dis_K[4];
	Result_Camera_Dis_K[5] = camera_Dis_K[5];
	delete[] camera_Dis_K;

	double* camera_Dis_P = new double[2];
	return_size = fread(camera_Dis_P, sizeof(double), 2, fp);
	if (return_size != 2)
	{
		QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Successfully, but no data"), QMessageBox::Ok, QMessageBox::Ok);
		clear_all_data();
		update_show_view();
		in_reading_keypoint = false;
		update_points_view();
		update_table_view();
		ui.open_path_lineEdit->setText(tr(""));
		fclose(fp);
		return;
	}
	Result_Camera_Dis_P[0] = camera_Dis_P[0];
	Result_Camera_Dis_P[1] = camera_Dis_P[1];
	delete[] camera_Dis_P;


	double* camera_Dis_T = new double[4];
	return_size = fread(camera_Dis_T, sizeof(double), 4, fp);
	if (return_size != 4)
	{
		QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Successfully, but no data"), QMessageBox::Ok, QMessageBox::Ok);
		clear_all_data();
		update_show_view();
		in_reading_keypoint = false;
		update_points_view();
		update_table_view();
		ui.open_path_lineEdit->setText(tr(""));
		fclose(fp);
		return;
	}
	Result_Camera_Dis_T[0] = camera_Dis_T[0];
	Result_Camera_Dis_T[1] = camera_Dis_T[1];
	Result_Camera_Dis_T[2] = camera_Dis_T[2];
	Result_Camera_Dis_T[3] = camera_Dis_T[3];
	delete[] camera_Dis_T;



	int image_size;
	return_size = fread(&image_size, sizeof(int), 1, fp);
	if (return_size != 1)
	{
		QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Successfully, but no data"), QMessageBox::Ok, QMessageBox::Ok);
		clear_all_data();
		update_show_view();
		in_reading_keypoint = false;
		update_points_view();
		update_table_view();
		ui.open_path_lineEdit->setText(tr(""));
		fclose(fp);
		return;
	}

	if (image_size == 0 || feof(fp))
	{
		QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Successfully, but no data"), QMessageBox::Ok, QMessageBox::Ok);
		clear_all_data();
		update_show_view();
		in_reading_keypoint = false;
		update_points_view();
		update_table_view();
		ui.open_path_lineEdit->setText(tr(""));
		fclose(fp);
		return;
	}
	for (int ii = 0; ii < image_size; ii++)
	{
		int image_index;
		return_size = fread(&image_index, sizeof(int), 1, fp);
		if (return_size != 1)
		{
			if (feof(fp))
			{
				break;
			}
			QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			clear_all_data();
			update_show_view();
			in_reading_keypoint = false;
			update_points_view();
			update_table_view();
			ui.open_path_lineEdit->setText(tr(""));
			fclose(fp);
			return;
		}
		long long string_path_length;
		return_size = fread(&string_path_length, sizeof(long long), 1, fp);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			clear_all_data();
			update_show_view();
			in_reading_keypoint = false;
			update_points_view();
			update_table_view();
			ui.open_path_lineEdit->setText(tr(""));
			fclose(fp);
			return;
		}
		std::string string_path;
		string_path.resize(string_path_length);
		return_size = fread(&string_path[0], sizeof(char), string_path_length, fp);
		if (return_size != string_path_length)
		{
			QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			clear_all_data();
			update_show_view();
			in_reading_keypoint = false;
			update_points_view();
			update_table_view();
			ui.open_path_lineEdit->setText(tr(""));
			fclose(fp);
			return;
		}
		Image_serial_name.push_back(QString(QString::fromLocal8Bit(string_path.data())));
		KeyPoint_Enable_serial.push_back(true);
		QDir file_temp = QDir(QString(QString::fromLocal8Bit(string_path.data())));
		file_temp.cdUp();
		ui.open_path_lineEdit->setText(file_temp.absolutePath());
		int image_width = -1, image_height = -1;
		return_size = fread(&image_width, sizeof(int), 1, fp);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			clear_all_data();
			update_show_view();
			in_reading_keypoint = false;
			update_points_view();
			update_table_view();
			ui.open_path_lineEdit->setText(tr(""));
			fclose(fp);
			return;
		}
		return_size = fread(&image_height, sizeof(int), 1, fp);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			clear_all_data();
			update_show_view();
			in_reading_keypoint = false;
			update_points_view();
			update_table_view();
			ui.open_path_lineEdit->setText(tr(""));
			fclose(fp);
			return;
		}
		ImageWidth_serial.push_back(image_width);
		ImageHeight_serial.push_back(image_height);


		double* camera_R_overall_temp = new double[4];
		return_size = fread(camera_R_overall_temp, sizeof(double), 4, fp);
		if (return_size != 4)
		{
			QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			clear_all_data();
			update_show_view();
			in_reading_keypoint = false;
			update_points_view();
			update_table_view();
			ui.open_path_lineEdit->setText(tr(""));
			fclose(fp);
			return;
		}
		Eigen::Quaterniond camera_R_overall(camera_R_overall_temp[0], camera_R_overall_temp[1], camera_R_overall_temp[2], camera_R_overall_temp[3]);
		delete[] camera_R_overall_temp;
		Result_R.push_back(camera_R_overall);

		double* camera_T_overall_temp = new double[3];
		return_size = fread(camera_T_overall_temp, sizeof(double), 3, fp);
		if (return_size != 3)
		{
			QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			clear_all_data();
			update_show_view();
			in_reading_keypoint = false;
			update_points_view();
			update_table_view();
			ui.open_path_lineEdit->setText(tr(""));
			fclose(fp);
			return;
		}
		Eigen::Vector3d camera_T_overall(camera_T_overall_temp[0], camera_T_overall_temp[1], camera_T_overall_temp[2]);
		delete[] camera_T_overall_temp;
		Result_T.push_back(camera_T_overall);

		int corn_size = -1;
		return_size = fread(&corn_size, sizeof(int), 1, fp);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			clear_all_data();
			update_show_view();
			in_reading_keypoint = false;
			update_points_view();
			update_table_view();
			ui.open_path_lineEdit->setText(tr(""));
			fclose(fp);
			return;
		}
		std::vector<Coded_detect_inf> kp_temp;
		std::vector<std::vector<cv::Point>> contours_temp;
		std::vector<cv::Point2d> error_now;
		for (int pp = 0; pp < corn_size; pp++)
		{
			double temp_x, temp_y, temp_ra, temp_rb, temp_theta, temp_err;
			double temp_N1, temp_N2, temp_N3, temp_N4, temp_N5, temp_N6;
			double reproj_x, reproj_y;
			int temp_code, temp_contours_num;
			return_size = fread(&temp_x, sizeof(double), 1, fp);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				clear_all_data();
				update_show_view();
				in_reading_keypoint = false;
				update_points_view();
				update_table_view();
				ui.open_path_lineEdit->setText(tr(""));
				fclose(fp);
				return;
			}
			return_size = fread(&temp_y, sizeof(double), 1, fp);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				clear_all_data();
				update_show_view();
				in_reading_keypoint = false;
				update_points_view();
				update_table_view();
				ui.open_path_lineEdit->setText(tr(""));
				fclose(fp);
				return;
			}
			return_size = fread(&temp_ra, sizeof(double), 1, fp);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				clear_all_data();
				update_show_view();
				in_reading_keypoint = false;
				update_points_view();
				update_table_view();
				ui.open_path_lineEdit->setText(tr(""));
				fclose(fp);
				return;
			}
			return_size = fread(&temp_rb, sizeof(double), 1, fp);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				clear_all_data();
				update_show_view();
				in_reading_keypoint = false;
				update_points_view();
				update_table_view();
				ui.open_path_lineEdit->setText(tr(""));
				fclose(fp);
				return;
			}
			return_size = fread(&temp_theta, sizeof(double), 1, fp);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				clear_all_data();
				update_show_view();
				in_reading_keypoint = false;
				update_points_view();
				update_table_view();
				ui.open_path_lineEdit->setText(tr(""));
				fclose(fp);
				return;
			}
			return_size = fread(&temp_code, sizeof(int), 1, fp);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				clear_all_data();
				update_show_view();
				in_reading_keypoint = false;
				update_points_view();
				update_table_view();
				ui.open_path_lineEdit->setText(tr(""));
				fclose(fp);
				return;
			}
			return_size = fread(&temp_err, sizeof(double), 1, fp);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				clear_all_data();
				update_show_view();
				in_reading_keypoint = false;
				update_points_view();
				update_table_view();
				ui.open_path_lineEdit->setText(tr(""));
				fclose(fp);
				return;
			}
			return_size = fread(&temp_N1, sizeof(double), 1, fp);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Read KeyPoint File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				clear_all_data();
				update_show_view();
				in_reading_keypoint = false;
				update_table_view();
				ui.open_path_lineEdit->setText(tr(""));
				fclose(fp);
				return;
			}
			return_size = fread(&temp_N2, sizeof(double), 1, fp);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Read KeyPoint File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				clear_all_data();
				update_show_view();
				in_reading_keypoint = false;
				update_table_view();
				ui.open_path_lineEdit->setText(tr(""));
				fclose(fp);
				return;
			}
			return_size = fread(&temp_N3, sizeof(double), 1, fp);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Read KeyPoint File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				clear_all_data();
				update_show_view();
				in_reading_keypoint = false;
				update_table_view();
				ui.open_path_lineEdit->setText(tr(""));
				fclose(fp);
				return;
			}
			return_size = fread(&temp_N4, sizeof(double), 1, fp);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Read KeyPoint File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				clear_all_data();
				update_show_view();
				in_reading_keypoint = false;
				update_table_view();
				ui.open_path_lineEdit->setText(tr(""));
				fclose(fp);
				return;
			}
			return_size = fread(&temp_N5, sizeof(double), 1, fp);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Read KeyPoint File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				clear_all_data();
				update_show_view();
				in_reading_keypoint = false;
				update_table_view();
				ui.open_path_lineEdit->setText(tr(""));
				fclose(fp);
				return;
			}
			return_size = fread(&temp_N6, sizeof(double), 1, fp);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Read KeyPoint File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				clear_all_data();
				update_show_view();
				in_reading_keypoint = false;
				update_table_view();
				ui.open_path_lineEdit->setText(tr(""));
				fclose(fp);
				return;
			}
			ellipse_view_pose temp_pose;
			temp_pose.view_N_1[0] = temp_N1;
			temp_pose.view_N_1[1] = temp_N2;
			temp_pose.view_N_1[2] = temp_N3;
			temp_pose.view_N_2[0] = temp_N4;
			temp_pose.view_N_2[1] = temp_N5;
			temp_pose.view_N_2[2] = temp_N6;
			Coded_detect_inf points_now_temp;
			points_now_temp.x = temp_x;
			points_now_temp.y = temp_y;
			points_now_temp.r_a = temp_ra;
			points_now_temp.r_b = temp_rb;
			points_now_temp.angle_in_pi = temp_theta;
			points_now_temp.code_num = temp_code;
			points_now_temp.fit_err = temp_err;
			points_now_temp.pose = temp_pose;
			kp_temp.push_back(points_now_temp);
			std::vector<cv::Point> contours_now;
			return_size = fread(&temp_contours_num, sizeof(int), 1, fp);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				clear_all_data();
				update_show_view();
				in_reading_keypoint = false;
				update_points_view();
				update_table_view();
				ui.open_path_lineEdit->setText(tr(""));
				fclose(fp);
				return;
			}
			for (int qq = 0; qq < temp_contours_num; qq++)
			{
				int temp_x_contour, temp_y_contour;
				return_size = fread(&temp_x_contour, sizeof(int), 1, fp);
				if (return_size != 1)
				{
					QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
					clear_all_data();
					update_show_view();
					in_reading_keypoint = false;
					update_points_view();
					update_table_view();
					ui.open_path_lineEdit->setText(tr(""));
					fclose(fp);
					return;
				}
				return_size = fread(&temp_y_contour, sizeof(int), 1, fp);
				if (return_size != 1)
				{
					QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
					clear_all_data();
					update_show_view();
					in_reading_keypoint = false;
					update_points_view();
					update_table_view();
					ui.open_path_lineEdit->setText(tr(""));
					fclose(fp);
					return;
				}
				contours_now.push_back(cv::Point(temp_x_contour, temp_y_contour));
			}
			contours_temp.push_back(contours_now);

			return_size = fread(&reproj_x, sizeof(double), 1, fp);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				clear_all_data();
				update_show_view();
				in_reading_keypoint = false;
				update_points_view();
				update_table_view();
				ui.open_path_lineEdit->setText(tr(""));
				fclose(fp);
				return;
			}

			return_size = fread(&reproj_y, sizeof(double), 1, fp);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				clear_all_data();
				update_show_view();
				in_reading_keypoint = false;
				update_points_view();
				update_table_view();
				ui.open_path_lineEdit->setText(tr(""));
				fclose(fp);
				return;
			}
			error_now.push_back(cv::Point2d(reproj_x, reproj_y));

		}
		Result_err_every.push_back(error_now);
		Result_contours_serial.push_back(contours_temp);
		KeyPoint_code_serial.push_back(kp_temp);
	}


	int points_size;
	return_size = fread(&points_size, sizeof(int), 1, fp);
	if (return_size != 1)
	{
		QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Successfully, but 3d points"), QMessageBox::Ok, QMessageBox::Ok);
		clear_all_data();
		update_show_view();
		in_reading_keypoint = false;
		update_points_view();
		update_table_view();
		ui.open_path_lineEdit->setText(tr(""));
		fclose(fp);
		return;
	}

	for (int ii = 0; ii < points_size; ii++)
	{
		int temp_code,temp_weight;
		double temp_x, temp_y, temp_z, temp_R;
		double temp_orix, temp_oriy, temp_oriz;
		return_size = fread(&temp_code, sizeof(int), 1, fp);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			clear_all_data();
			update_show_view();
			in_reading_keypoint = false;
			update_points_view();
			update_table_view();
			ui.open_path_lineEdit->setText(tr(""));
			fclose(fp);
			return;
		}
		return_size = fread(&temp_weight, sizeof(int), 1, fp);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			clear_all_data();
			update_show_view();
			in_reading_keypoint = false;
			update_points_view();
			update_table_view();
			ui.open_path_lineEdit->setText(tr(""));
			fclose(fp);
			return;
		}
		return_size = fread(&temp_x, sizeof(double), 1, fp);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			clear_all_data();
			update_show_view();
			in_reading_keypoint = false;
			update_points_view();
			update_table_view();
			ui.open_path_lineEdit->setText(tr(""));
			fclose(fp);
			return;
		}
		return_size = fread(&temp_y, sizeof(double), 1, fp);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			clear_all_data();
			update_show_view();
			in_reading_keypoint = false;
			update_points_view();
			update_table_view();
			ui.open_path_lineEdit->setText(tr(""));
			fclose(fp);
			return;
		}
		return_size = fread(&temp_z, sizeof(double), 1, fp);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			clear_all_data();
			update_show_view();
			in_reading_keypoint = false;
			update_points_view();
			update_table_view();
			ui.open_path_lineEdit->setText(tr(""));
			fclose(fp);
			return;
		}
		return_size = fread(&temp_R, sizeof(double), 1, fp);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			clear_all_data();
			update_show_view();
			in_reading_keypoint = false;
			update_points_view();
			update_table_view();
			ui.open_path_lineEdit->setText(tr(""));
			fclose(fp);
			return;
		}
		return_size = fread(&temp_orix, sizeof(double), 1, fp);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			clear_all_data();
			update_show_view();
			in_reading_keypoint = false;
			update_points_view();
			update_table_view();
			ui.open_path_lineEdit->setText(tr(""));
			fclose(fp);
			return;
		}
		return_size = fread(&temp_oriy, sizeof(double), 1, fp);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			clear_all_data();
			update_show_view();
			in_reading_keypoint = false;
			update_points_view();
			update_table_view();
			ui.open_path_lineEdit->setText(tr(""));
			fclose(fp);
			return;
		}
		return_size = fread(&temp_oriz, sizeof(double), 1, fp);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			clear_all_data();
			update_show_view();
			in_reading_keypoint = false;
			update_points_view();
			update_table_view();
			ui.open_path_lineEdit->setText(tr(""));
			fclose(fp);
			return;
		}

		Result_xyz_serial.push_back(sfm_3d_group(temp_code, temp_weight, temp_x, temp_y, temp_z, 0, temp_R));
		xyz_need_levelling_serial.push_back(false);
		xyz_need_levelling_code_serial.push_back(temp_code);
		Result_Ori.push_back(Eigen::Vector3d(temp_orix, temp_oriy, temp_oriz));
	}
	std::sort(xyz_need_levelling_code_serial.begin(), xyz_need_levelling_code_serial.end());

	int queue_size;
	return_size = fread(&queue_size, sizeof(int), 1, fp);
	if (return_size != 1)
	{
		QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Successfully, but 3d points"), QMessageBox::Ok, QMessageBox::Ok);
		clear_all_data();
		update_show_view();
		in_reading_keypoint = false;
		update_points_view();
		update_table_view();
		ui.open_path_lineEdit->setText(tr(""));
		fclose(fp);
		return;
	}

	for (int ii = 0; ii < queue_size; ii++)
	{
		int temp_code;
		return_size = fread(&temp_code, sizeof(int), 1, fp);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			clear_all_data();
			update_show_view();
			in_reading_keypoint = false;
			update_points_view();
			update_table_view();
			ui.open_path_lineEdit->setText(tr(""));
			fclose(fp);
			return;
		}
		Keypoints_index_in_queue.push_back(temp_code);
	}

	int scale_size;
	return_size = fread(&scale_size, sizeof(int), 1, fp);
	if (return_size != 1)
	{
		QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Successfully, but 3d points"), QMessageBox::Ok, QMessageBox::Ok);
		clear_all_data();
		update_show_view();
		in_reading_keypoint = false;
		update_points_view();
		update_table_view();
		ui.open_path_lineEdit->setText(tr(""));
		fclose(fp);
		return;
	}

	for (int ii = 0; ii < scale_size; ii++)
	{
		add_scale_lable();
		int begin_code,end_code;
		double dis_value;
		return_size = fread(&begin_code, sizeof(int), 1, fp);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			clear_all_data();
			update_show_view();
			in_reading_keypoint = false;
			update_points_view();
			update_table_view();
			ui.open_path_lineEdit->setText(tr(""));
			fclose(fp);
			return;
		}
		return_size = fread(&end_code, sizeof(int), 1, fp);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			clear_all_data();
			update_show_view();
			in_reading_keypoint = false;
			update_points_view();
			update_table_view();
			ui.open_path_lineEdit->setText(tr(""));
			fclose(fp);
			return;
		}
		return_size = fread(&dis_value, sizeof(double), 1, fp);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			clear_all_data();
			update_show_view();
			in_reading_keypoint = false;
			update_points_view();
			update_table_view();
			ui.open_path_lineEdit->setText(tr(""));
			fclose(fp);
			return;
		}
		begin_scale_box[begin_scale_box.size() - 1]->setValue(begin_code);
		end_scale_box[end_scale_box.size() - 1]->setValue(end_code);
		value_scale_box[value_scale_box.size() - 1]->setValue(dis_value);
	}

	int fit_size;
	return_size = fread(&fit_size, sizeof(int), 1, fp);
	if (return_size != 1)
	{
		QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Successfully, but 3d points"), QMessageBox::Ok, QMessageBox::Ok);
		clear_all_data();
		update_show_view();
		in_reading_keypoint = false;
		update_points_view();
		update_table_view();
		ui.open_path_lineEdit->setText(tr(""));
		fclose(fp);
		return;
	}

	for (int ii = 0; ii < fit_size; ii++)
	{
		int begin_code;
		return_size = fread(&begin_code, sizeof(int), 1, fp);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Read SFM Result File"), tr("Read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			clear_all_data();
			update_show_view();
			in_reading_keypoint = false;
			update_points_view();
			update_table_view();
			ui.open_path_lineEdit->setText(tr(""));
			fclose(fp);
			return;
		}
		for (int jj = 0; jj < xyz_need_levelling_code_serial.size(); jj++)
		{
			if (xyz_need_levelling_code_serial[jj] == begin_code)
			{
				xyz_need_levelling_serial[jj] = true;
			}
		}
	}
	for (int ii = 0; ii < Keypoints_index_in_queue.size(); ii++)
	{
		if (Keypoints_index_in_queue[ii] >= Result_R.size() || Keypoints_index_in_queue[ii] >= Result_T.size())
		{
			QMessageBox::warning(this, tr("Read SFM Result File"), tr("Calculation sequence logic error, read Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			clear_all_data();
			update_show_view();
			in_reading_keypoint = false;
			update_points_view();
			update_table_view();
			ui.open_path_lineEdit->setText(tr(""));
			fclose(fp);
			return;
		}
		cv::Mat mat_temp;
		cv::eigen2cv(Result_R[ii].matrix(), mat_temp);
		cv::Mat mat_temp2;
		cv::eigen2cv(Result_T[ii], mat_temp2);
		Keypoints_R_in_queue.push_back(mat_temp);
		Keypoints_T_in_queue.push_back(mat_temp2);
	}
	if (Image_serial_name.size() == 0)
	{
		ui.show_detect_spinBox->setEnabled(false);
		ui.show_detect_horizontalSlider->setEnabled(false);
		ui.show_detect_spinBox->setMinimum(1);
		ui.show_detect_spinBox->setMaximum(1);
		ui.show_detect_horizontalSlider->setMinimum(1);
		ui.show_detect_horizontalSlider->setMaximum(1);
		ui.show_detect_spinBox->setValue(1);
		ui.show_detect_horizontalSlider->setValue(1);
		ui.Detect_view->clear_image();
	}
	else
	{
		ui.show_detect_spinBox->setEnabled(true);
		ui.show_detect_horizontalSlider->setEnabled(true);
		ui.show_detect_spinBox->setMinimum(1);
		ui.show_detect_spinBox->setMaximum(Image_serial_name.size());
		ui.show_detect_horizontalSlider->setMinimum(1);
		ui.show_detect_horizontalSlider->setMaximum(Image_serial_name.size());
		ui.show_detect_spinBox->setValue(1);
		ui.show_detect_horizontalSlider->setValue(1);
	}

	in_reading_keypoint = false;
	need_recal_keypoints = false;
	update_show_view();
	update_points_view();
	update_table_view();
	QMessageBox::information(this, tr("Read SFM Result File"), tr("Read Successfully."), QMessageBox::Ok);
	fclose(fp);
}

void SFM_module::save_KeyPoint_file(QString path_save)
{
	FILE* f = fopen(path_save.toLocal8Bit().data(), "wb");
	unsigned char save_mode = 31;
	size_t return_size = fwrite(&save_mode, sizeof(unsigned char), 1, f);
	if (return_size != 1)
	{
		QMessageBox::warning(this, tr("Save KeyPoint File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
		fclose(f);
		return;
	}
	for (int jj = 0; jj < Image_serial_name.size(); jj++)
	{
		int image_index = jj;
		return_size = fwrite(&image_index, sizeof(int), 1, f);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Save KeyPoint File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			fclose(f);
			return;
		}
		long long string_len = Image_serial_name[jj].toLocal8Bit().toStdString().size();
		return_size = fwrite(&string_len, sizeof(long long), 1, f);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Save KeyPoint File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			fclose(f);
			return;
		}
		std::string im_path = Image_serial_name[jj].toLocal8Bit().toStdString();
		return_size = fwrite(&im_path[0], sizeof(char), string_len, f);
		if (return_size != string_len)
		{
			QMessageBox::warning(this, tr("Save KeyPoint File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			fclose(f);
			return;
		}
		return_size = fwrite(&ImageWidth_serial[jj], sizeof(int), 1, f);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Save KeyPoint File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			fclose(f);
			return;
		}
		return_size = fwrite(&ImageHeight_serial[jj], sizeof(int), 1, f);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Save KeyPoint File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			fclose(f);
			return;
		}
		int corn_number = KeyPoint_code_serial[jj].size();
		return_size = fwrite(&corn_number, sizeof(int), 1, f);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Save KeyPoint File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			fclose(f);
			return;
		}
		for (int pp = 0; pp < corn_number; pp++)
		{
			return_size = fwrite(&(KeyPoint_code_serial_update[jj][pp].x), sizeof(double), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save KeyPoint File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
			return_size = fwrite(&(KeyPoint_code_serial_update[jj][pp].y), sizeof(double), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save KeyPoint File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
			return_size = fwrite(&(KeyPoint_code_serial_update[jj][pp].r_a), sizeof(double), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save KeyPoint File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
			return_size = fwrite(&(KeyPoint_code_serial_update[jj][pp].r_b), sizeof(double), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save KeyPoint File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
			return_size = fwrite(&(KeyPoint_code_serial_update[jj][pp].angle_in_pi), sizeof(double), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save KeyPoint File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
			return_size = fwrite(&(KeyPoint_code_serial_update[jj][pp].code_num), sizeof(int), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save KeyPoint File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
			return_size = fwrite(&(KeyPoint_code_serial_update[jj][pp].fit_err), sizeof(double), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save KeyPoint File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
			return_size = fwrite(&(KeyPoint_code_serial_update[jj][pp].pose.view_N_1[0]), sizeof(double), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save KeyPoint File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
			return_size = fwrite(&(KeyPoint_code_serial_update[jj][pp].pose.view_N_1[1]), sizeof(double), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save KeyPoint File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
			return_size = fwrite(&(KeyPoint_code_serial_update[jj][pp].pose.view_N_1[2]), sizeof(double), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save KeyPoint File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
			return_size = fwrite(&(KeyPoint_code_serial_update[jj][pp].pose.view_N_2[0]), sizeof(double), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save KeyPoint File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
			return_size = fwrite(&(KeyPoint_code_serial_update[jj][pp].pose.view_N_2[1]), sizeof(double), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save KeyPoint File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
			return_size = fwrite(&(KeyPoint_code_serial_update[jj][pp].pose.view_N_2[2]), sizeof(double), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save KeyPoint File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
			int contours_number = Result_contours_serial[jj][pp].size();
			return_size = fwrite(&(contours_number), sizeof(int), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save KeyPoint File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
			for (int qq = 0; qq < contours_number; qq++)
			{
				return_size = fwrite(&(Result_contours_serial[jj][pp][qq].x), sizeof(int), 1, f);
				if (return_size != 1)
				{
					QMessageBox::warning(this, tr("Save KeyPoint File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
					fclose(f);
					return;
				}
				return_size = fwrite(&(Result_contours_serial[jj][pp][qq].y), sizeof(int), 1, f);
				if (return_size != 1)
				{
					QMessageBox::warning(this, tr("Save KeyPoint File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
					fclose(f);
					return;
				}
			}
		}
	}
	QMessageBox::information(this, tr("Save KeyPoint File"), tr("Save Successfully."), QMessageBox::Ok);
	fclose(f);
}
void SFM_module::save_KeyPoint_file_csv(QString path_save)
{
	QFile file(path_save);
	if (!file.exists())  //文件不存在的时新建
	{
		file.open(QIODevice::WriteOnly | QIODevice::Text);
		QTextStream txtOutPut(&file);
		file.close();
	}
	file.open(QIODevice::WriteOnly | QIODevice::Truncate | QIODevice::Text);
	QTextStream txtOutPut(&file);
	txtOutPut.setRealNumberPrecision(16);
	txtOutPut.setEncoding(QStringConverter::System);
	for (int jj = 0; jj < Image_serial_name.size(); jj++)
	{
		txtOutPut <<tr("Camera-Index:") << jj + 1 << ",";
		txtOutPut << Image_serial_name[jj].toLocal8Bit() << ",";
		txtOutPut << ImageWidth_serial[jj] << ",";
		txtOutPut << ImageHeight_serial[jj];

		txtOutPut << "\n" << tr("Code");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "," << KeyPoint_code_serial_update[jj][pp].code_num;
		}
		txtOutPut << "\n" << tr("X_pixel");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "," << KeyPoint_code_serial_update[jj][pp].x;
		}
		txtOutPut << "\n" << tr("Y_pixel");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "," << KeyPoint_code_serial_update[jj][pp].y;
		}
		txtOutPut << "\n" << tr("R_a");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "," << KeyPoint_code_serial_update[jj][pp].r_a;
		}
		txtOutPut << "\n" << tr("R_b");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "," << KeyPoint_code_serial_update[jj][pp].r_b;
		}
		txtOutPut << "\n" << tr("Angle");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "," << KeyPoint_code_serial_update[jj][pp].angle_in_pi;
		}
		txtOutPut << "\n" << tr("Pose_1_Nx");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "," << KeyPoint_code_serial_update[jj][pp].pose.view_N_1[0];
		}
		txtOutPut << "\n" << tr("Pose_1_Ny");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "," << KeyPoint_code_serial_update[jj][pp].pose.view_N_1[1];
		}
		txtOutPut << "\n" << tr("Pose_1_Nz");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "," << KeyPoint_code_serial_update[jj][pp].pose.view_N_1[2];
		}
		txtOutPut << "\n" << tr("Pose_2_Nx");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "," << KeyPoint_code_serial_update[jj][pp].pose.view_N_2[0];
		}
		txtOutPut << "\n" << tr("Pose_2_Ny");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "," << KeyPoint_code_serial_update[jj][pp].pose.view_N_2[1];
		}
		txtOutPut << "\n" << tr("Pose_2_Nz");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "," << KeyPoint_code_serial_update[jj][pp].pose.view_N_2[2];
		}
		txtOutPut << "\n";
	}
	//sbsbsb
	QMessageBox::information(this, tr("Save KeyPoint File"), tr("Save Successfully."), QMessageBox::Ok);
	file.close();
}
void SFM_module::save_KeyPoint_file_txt(QString path_save)
{
	QFile file(path_save);
	if (!file.exists())  //文件不存在的时新建
	{
		file.open(QIODevice::WriteOnly | QIODevice::Text);
		QTextStream txtOutPut(&file);
		file.close();
	}
	file.open(QIODevice::WriteOnly | QIODevice::Truncate | QIODevice::Text);
	QTextStream txtOutPut(&file);
	txtOutPut.setRealNumberPrecision(16);
	txtOutPut.setEncoding(QStringConverter::System);
	for (int jj = 0; jj < Image_serial_name.size(); jj++)
	{
		txtOutPut << tr("Camera-Index:") << jj + 1 << "\t";
		txtOutPut << Image_serial_name[jj].toLocal8Bit() << "\t";
		txtOutPut << ImageWidth_serial[jj] << "\t";
		txtOutPut << ImageHeight_serial[jj];

		txtOutPut << "\n" << tr("Code");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "\t" << KeyPoint_code_serial_update[jj][pp].code_num;
		}
		txtOutPut << "\n" << tr("X_pixel");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "\t" << KeyPoint_code_serial_update[jj][pp].x;
		}
		txtOutPut << "\n" << tr("Y_pixel");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "\t" << KeyPoint_code_serial_update[jj][pp].y;
		}
		txtOutPut << "\n" << tr("R_a");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "\t" << KeyPoint_code_serial_update[jj][pp].r_a;
		}
		txtOutPut << "\n" << tr("R_b");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "\t" << KeyPoint_code_serial_update[jj][pp].r_b;
		}
		txtOutPut << "\n" << tr("Angle");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "\t" << KeyPoint_code_serial_update[jj][pp].angle_in_pi;
		}
		txtOutPut << "\n" << tr("Pose_1_Nx");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "\t" << KeyPoint_code_serial_update[jj][pp].pose.view_N_1[0];
		}
		txtOutPut << "\n" << tr("Pose_1_Ny");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "\t" << KeyPoint_code_serial_update[jj][pp].pose.view_N_1[1];
		}
		txtOutPut << "\n" << tr("Pose_1_Nz");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "\t" << KeyPoint_code_serial_update[jj][pp].pose.view_N_1[2];
		}
		txtOutPut << "\n" << tr("Pose_2_Nx");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "\t" << KeyPoint_code_serial_update[jj][pp].pose.view_N_2[0];
		}
		txtOutPut << "\n" << tr("Pose_2_Ny");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "\t" << KeyPoint_code_serial_update[jj][pp].pose.view_N_2[1];
		}
		txtOutPut << "\n" << tr("Pose_2_Nz");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "\t" << KeyPoint_code_serial_update[jj][pp].pose.view_N_2[2];
		}
		txtOutPut << "\n";
	}
	QMessageBox::information(this, tr("Save KeyPoint File"), tr("Save Successfully."), QMessageBox::Ok);
	file.close();
}
void SFM_module::save_KeyPoint_file_mat(QString path_save)
{
	MATFile* pmatFile = NULL;
	pmatFile = matOpen(path_save.toLocal8Bit().data(), "w");
	if (pmatFile == NULL)
	{
		QMessageBox::warning(this, tr("Save KeyPoint File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
		return;
	}

	for (int jj = 0; jj < Image_serial_name.size(); jj++)
	{
		mxArray* pMxArray = NULL;
		pMxArray = mxCreateDoubleMatrix(12, KeyPoint_code_serial[jj].size(), mxREAL);
		if (pMxArray == NULL)
		{
			QMessageBox::warning(this, tr("Save KeyPoint File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			return;
		}
		double* out_data = new double[KeyPoint_code_serial[jj].size() * 12];
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			out_data[pp * 12] = KeyPoint_code_serial[jj][pp].code_num;
			out_data[pp * 12 + 1] = KeyPoint_code_serial[jj][pp].x;
			out_data[pp * 12 + 2] = KeyPoint_code_serial[jj][pp].y;
			out_data[pp * 12 + 3] = KeyPoint_code_serial[jj][pp].r_a;
			out_data[pp * 12 + 4] = KeyPoint_code_serial[jj][pp].r_b;
			out_data[pp * 12 + 5] = KeyPoint_code_serial[jj][pp].angle_in_pi;
			out_data[pp * 12 + 6] = KeyPoint_code_serial[jj][pp].pose.view_N_1[0];
			out_data[pp * 12 + 7] = KeyPoint_code_serial[jj][pp].pose.view_N_1[1];
			out_data[pp * 12 + 8] = KeyPoint_code_serial[jj][pp].pose.view_N_1[2];
			out_data[pp * 12 + 9] = KeyPoint_code_serial[jj][pp].pose.view_N_2[0];
			out_data[pp * 12 + 10] = KeyPoint_code_serial[jj][pp].pose.view_N_2[1];
			out_data[pp * 12 + 11] = KeyPoint_code_serial[jj][pp].pose.view_N_2[2];
		}
		memcpy((void*)(mxGetPr(pMxArray)), (void*)out_data, sizeof(double)* KeyPoint_code_serial[jj].size() * 12);
		QString var_mat_name = tr("Image_index") + QString("%1").arg(jj + 1, 8, 10, QLatin1Char('0'));
		matPutVariable(pmatFile, var_mat_name.toLocal8Bit().toStdString().data(), pMxArray);
		delete[] out_data;
		mxDestroyArray(pMxArray);
	}
	QMessageBox::information(this, tr("Save KeyPoint File"), tr("Save Successfully."), QMessageBox::Ok);
	matClose(pmatFile);
}
void SFM_module::save_SFMResult_file(QString path_save)
{
	FILE* f = fopen(path_save.toLocal8Bit().data(), "wb");
	unsigned char save_mode = 32;
	size_t return_size = fwrite(&save_mode, sizeof(unsigned char), 1, f);
	if (return_size != 1)
	{
		QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
		fclose(f);
		return;
	}

	return_size = fwrite(Result_Camera_K, sizeof(double), 5, f);
	if (return_size != 5)
	{
		QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
		fclose(f);
		return;
	}

	return_size = fwrite(Result_Camera_Dis_K, sizeof(double), 6, f);
	if (return_size != 6)
	{
		QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
		fclose(f);
		return;
	}

	return_size = fwrite(Result_Camera_Dis_P, sizeof(double), 2, f);
	if (return_size != 2)
	{
		QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
		fclose(f);
		return;
	}

	return_size = fwrite(Result_Camera_Dis_T, sizeof(double), 4, f);
	if (return_size != 4)
	{
		QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
		fclose(f);
		return;
	}

	int image_size = Image_serial_name.size();
	return_size = fwrite(&image_size, sizeof(int), 1, f);
	if (return_size != 1)
	{
		QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
		fclose(f);
		return;
	}

	for (int jj = 0; jj < Image_serial_name.size(); jj++)
	{
		int image_index = jj;
		return_size = fwrite(&image_index, sizeof(int), 1, f);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			fclose(f);
			return;
		}
		long long string_len = Image_serial_name[jj].toLocal8Bit().toStdString().size();
		return_size = fwrite(&string_len, sizeof(long long), 1, f);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			fclose(f);
			return;
		}
		std::string im_path = Image_serial_name[jj].toLocal8Bit().toStdString();
		return_size = fwrite(&im_path[0], sizeof(char), string_len, f);
		if (return_size != string_len)
		{
			QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			fclose(f);
			return;
		}
		return_size = fwrite(&ImageWidth_serial[jj], sizeof(int), 1, f);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			fclose(f);
			return;
		}
		return_size = fwrite(&ImageHeight_serial[jj], sizeof(int), 1, f);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			fclose(f);
			return;
		}
		double* R_temp = new double[4];
		R_temp[0] = Result_R[jj].x();
		R_temp[1] = Result_R[jj].y();
		R_temp[2] = Result_R[jj].z();
		R_temp[3] = Result_R[jj].w();
		return_size = fwrite(R_temp, sizeof(double), 4, f);
		if (return_size != 4)
		{
			QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			fclose(f);
			return;
		}
		delete[] R_temp;

		double* T_temp = new double[3];
		T_temp[0] = Result_T[jj].x();
		T_temp[1] = Result_T[jj].y();
		T_temp[2] = Result_T[jj].z();
		return_size = fwrite(T_temp, sizeof(double), 3, f);
		if (return_size != 3)
		{
			QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			fclose(f);
			return;
		}
		delete[] T_temp;
		int corn_number = KeyPoint_code_serial[jj].size();
		return_size = fwrite(&corn_number, sizeof(int), 1, f);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			fclose(f);
			return;
		}
		for (int pp = 0; pp < corn_number; pp++)
		{
			return_size = fwrite(&(KeyPoint_code_serial_update[jj][pp].x), sizeof(double), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
			return_size = fwrite(&(KeyPoint_code_serial_update[jj][pp].y), sizeof(double), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
			return_size = fwrite(&(KeyPoint_code_serial_update[jj][pp].r_a), sizeof(double), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
			return_size = fwrite(&(KeyPoint_code_serial_update[jj][pp].r_b), sizeof(double), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
			return_size = fwrite(&(KeyPoint_code_serial_update[jj][pp].angle_in_pi), sizeof(double), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
			return_size = fwrite(&(KeyPoint_code_serial_update[jj][pp].code_num), sizeof(int), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
			return_size = fwrite(&(KeyPoint_code_serial_update[jj][pp].fit_err), sizeof(double), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
			return_size = fwrite(&(KeyPoint_code_serial_update[jj][pp].pose.view_N_1[0]), sizeof(double), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save KeyPoint File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
			return_size = fwrite(&(KeyPoint_code_serial_update[jj][pp].pose.view_N_1[1]), sizeof(double), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save KeyPoint File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
			return_size = fwrite(&(KeyPoint_code_serial_update[jj][pp].pose.view_N_1[2]), sizeof(double), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save KeyPoint File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
			return_size = fwrite(&(KeyPoint_code_serial_update[jj][pp].pose.view_N_2[0]), sizeof(double), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save KeyPoint File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
			return_size = fwrite(&(KeyPoint_code_serial_update[jj][pp].pose.view_N_2[1]), sizeof(double), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save KeyPoint File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
			return_size = fwrite(&(KeyPoint_code_serial_update[jj][pp].pose.view_N_2[2]), sizeof(double), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save KeyPoint File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
			int contours_number = Result_contours_serial[jj][pp].size();
			return_size = fwrite(&(contours_number), sizeof(int), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
			for (int qq = 0; qq < contours_number; qq++)
			{
				return_size = fwrite(&(Result_contours_serial[jj][pp][qq].x), sizeof(int), 1, f);
				if (return_size != 1)
				{
					QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
					fclose(f);
					return;
				}
				return_size = fwrite(&(Result_contours_serial[jj][pp][qq].y), sizeof(int), 1, f);
				if (return_size != 1)
				{
					QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
					fclose(f);
					return;
				}
			}
			return_size = fwrite(&(Result_err_every[jj][pp].x), sizeof(double), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
			return_size = fwrite(&(Result_err_every[jj][pp].y), sizeof(double), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
		}
	}

	int xyz_number = Result_xyz_serial.size();
	return_size = fwrite(&xyz_number, sizeof(int), 1, f);
	if (return_size != 1)
	{
		QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
		fclose(f);
		return;
	}
	for (int pp = 0; pp < xyz_number; pp++)
	{
		return_size = fwrite(&(Result_xyz_serial[pp].code), sizeof(int), 1, f);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			fclose(f);
			return;
		}
		return_size = fwrite(&(Result_xyz_serial[pp].weight), sizeof(int), 1, f);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			fclose(f);
			return;
		}
		return_size = fwrite(&(Result_xyz_serial[pp].x), sizeof(double), 1, f);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			fclose(f);
			return;
		}
		return_size = fwrite(&(Result_xyz_serial[pp].y), sizeof(double), 1, f);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			fclose(f);
			return;
		}
		return_size = fwrite(&(Result_xyz_serial[pp].z), sizeof(double), 1, f);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			fclose(f);
			return;
		}
		return_size = fwrite(&(Result_xyz_serial[pp].R), sizeof(double), 1, f);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			fclose(f);
			return;
		}
		if (Result_Ori.size() > pp)
		{
			return_size = fwrite(&(Result_Ori[pp].x()), sizeof(double), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
			return_size = fwrite(&(Result_Ori[pp].y()), sizeof(double), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
			return_size = fwrite(&(Result_Ori[pp].z()), sizeof(double), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
		}
		else
		{
			double ori_x = 0, ori_y = 0, ori_z = 0;
			return_size = fwrite(&(ori_x), sizeof(double), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
			return_size = fwrite(&(ori_y), sizeof(double), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
			return_size = fwrite(&(ori_z), sizeof(double), 1, f);
			if (return_size != 1)
			{
				QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
				fclose(f);
				return;
			}
		}
	}

	int index_queue_number = Keypoints_index_in_queue.size();
	return_size = fwrite(&index_queue_number, sizeof(int), 1, f);
	if (return_size != 1)
	{
		QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
		fclose(f);
		return;
	}
	for (int pp = 0; pp < index_queue_number; pp++)
	{
		return_size = fwrite(&(Keypoints_index_in_queue[pp]), sizeof(int), 1, f);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			fclose(f);
			return;
		}
	}

	int scale_number = begin_scale_box.size();
	return_size = fwrite(&scale_number, sizeof(int), 1, f);
	if (return_size != 1)
	{
		QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
		fclose(f);
		return;
	}
	for (int pp = 0; pp < scale_number; pp++)
	{
		int begin_code = begin_scale_box[pp]->value();
		int end_code = end_scale_box[pp]->value();
		double dis_value = value_scale_box[pp]->value();
		return_size = fwrite(&(begin_code), sizeof(int), 1, f);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			fclose(f);
			return;
		}
		return_size = fwrite(&(end_code), sizeof(int), 1, f);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			fclose(f);
			return;
		}
		return_size = fwrite(&(dis_value), sizeof(double), 1, f);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			fclose(f);
			return;
		}
	}

	int fit_number = 0;
	std::vector<int> fit_code;
	for (int ii = 0; ii < xyz_need_levelling_serial.size(); ii++)
	{
		if (xyz_need_levelling_serial[ii])
		{
			fit_code.push_back(xyz_need_levelling_code_serial[ii]);
			fit_number++;
		}
	}
	return_size = fwrite(&fit_number, sizeof(int), 1, f);
	if (return_size != 1)
	{
		QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
		fclose(f);
		return;
	}
	for (int pp = 0; pp < fit_number; pp++)
	{
		return_size = fwrite(&(fit_code[pp]), sizeof(int), 1, f);
		if (return_size != 1)
		{
			QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			fclose(f);
			return;
		}
	}

	QMessageBox::information(this, tr("Save Result File"), tr("Save Successfully."), QMessageBox::Ok);
	fclose(f);
}
void SFM_module::save_SFMResult_file_csv(QString path_save)
{
	QFile file(path_save);
	if (!file.exists())  //文件不存在的时新建
	{
		file.open(QIODevice::WriteOnly | QIODevice::Text);
		QTextStream txtOutPut(&file);
		file.close();
	}
	file.open(QIODevice::WriteOnly | QIODevice::Truncate | QIODevice::Text);
	QTextStream txtOutPut(&file);
	txtOutPut.setRealNumberPrecision(16);
	txtOutPut.setEncoding(QStringConverter::System);

	txtOutPut << tr("Camera Parameter Information") << "\n";
	txtOutPut << tr("IntrinsicMatrix") << "," << Result_Camera_K[0] << "," << Result_Camera_K[4] << "," << Result_Camera_K[2] << "\n";
	txtOutPut << tr("") << "," << 0 << "," << Result_Camera_K[1] << "," << Result_Camera_K[3]  << "\n";
	txtOutPut << tr("") << "," << 0 << "," << 0 << "," << 1 << "\n";
	txtOutPut << tr("RadialDistortion") << "," << Result_Camera_Dis_K[0] << "," << Result_Camera_Dis_K[1] << "," << Result_Camera_Dis_K[2] << ","
		<< Result_Camera_Dis_K[3] << "," << Result_Camera_Dis_K[4] << "," << Result_Camera_Dis_K[5]
		<< "\n";
	txtOutPut << tr("TangentialDistortion") << "," << Result_Camera_Dis_P[0] << "," << Result_Camera_Dis_P[1] << "\n";
	txtOutPut << tr("ThinPrismDistortion") << "," << Result_Camera_Dis_T[0] << "," << Result_Camera_Dis_T[1] 
		<< "," << Result_Camera_Dis_T[2] << "," << Result_Camera_Dis_T[3] << "\n";

	txtOutPut << "\n" << "\n" << "\n";

	txtOutPut << tr("3D points Information");
	txtOutPut << "\n" << tr("Code");
	for (int pp = 0; pp < Result_xyz_after_op_serial.size(); pp++)
	{
		txtOutPut << "," << Result_xyz_after_op_serial[pp].code;
	}
	txtOutPut << "\n" << tr("View Count");
	for (int pp = 0; pp < Result_xyz_after_op_serial.size(); pp++)
	{
		txtOutPut << "," << Result_xyz_after_op_serial[pp].weight;
	}
	txtOutPut << "\n" << tr("X");
	for (int pp = 0; pp < Result_xyz_after_op_serial.size(); pp++)
	{
		txtOutPut << "," << Result_xyz_after_op_serial[pp].x;
	}
	txtOutPut << "\n" << tr("Y");
	for (int pp = 0; pp < Result_xyz_after_op_serial.size(); pp++)
	{
		txtOutPut << "," << Result_xyz_after_op_serial[pp].y;
	}
	txtOutPut << "\n" << tr("Z");
	for (int pp = 0; pp < Result_xyz_after_op_serial.size(); pp++)
	{
		txtOutPut << "," << Result_xyz_after_op_serial[pp].z;
	}
	txtOutPut << "\n" << tr("R");
	for (int pp = 0; pp < Result_xyz_after_op_serial.size(); pp++)
	{
		txtOutPut << "," << Result_xyz_after_op_serial[pp].R;
	}
	txtOutPut << "\n" << tr("Ori_x");
	for (int pp = 0; pp < Result_xyz_after_op_serial.size(); pp++)
	{
		if (Result_Ori_after_op.size() > pp)
		{
			txtOutPut << "," << Result_Ori_after_op[pp].x();
		}
		else
		{
			txtOutPut << "," << tr("NAN");
		}
	}
	txtOutPut << "\n" << tr("Ori_y");
	for (int pp = 0; pp < Result_xyz_after_op_serial.size(); pp++)
	{
		if (Result_Ori_after_op.size() > pp)
		{
			txtOutPut << "," << Result_Ori_after_op[pp].y();
		}
		else
		{
			txtOutPut << "," << tr("NAN");
		}
	}
	txtOutPut << "\n" << tr("Ori_z");
	for (int pp = 0; pp < Result_xyz_after_op_serial.size(); pp++)
	{
		if (Result_Ori_after_op.size() > pp)
		{
			txtOutPut << "," << Result_Ori_after_op[pp].z();
		}
		else
		{
			txtOutPut << "," << tr("NAN");
		}
	}
	txtOutPut << "\n" << "\n" << "\n";

	txtOutPut << tr("Image-Keypoint Information") << "\n";

	for (int jj = 0; jj < Image_serial_name.size(); jj++)
	{
		txtOutPut << tr("Camera-Index:") << jj + 1 << ",";
		txtOutPut << Image_serial_name[jj].toLocal8Bit() << ",";
		txtOutPut << ImageWidth_serial[jj] << ",";
		txtOutPut << ImageHeight_serial[jj];

		txtOutPut << "\n" << tr("Code");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "," << KeyPoint_code_serial_update[jj][pp].code_num;
		}
		txtOutPut << "\n" << tr("X_pixel");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "," << KeyPoint_code_serial_update[jj][pp].x;
		}
		txtOutPut << "\n" << tr("Y_pixel");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "," << KeyPoint_code_serial_update[jj][pp].y;
		}
		txtOutPut << "\n" << tr("R_a");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "," << KeyPoint_code_serial_update[jj][pp].r_a;
		}
		txtOutPut << "\n" << tr("R_b");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "," << KeyPoint_code_serial_update[jj][pp].r_b;
		}
		txtOutPut << "\n" << tr("Angle");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "," << KeyPoint_code_serial_update[jj][pp].angle_in_pi;
		}
		txtOutPut << "\n" << tr("Pose_1_Nx");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "," << KeyPoint_code_serial_update[jj][pp].pose.view_N_1[0];
		}
		txtOutPut << "\n" << tr("Pose_1_Ny");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "," << KeyPoint_code_serial_update[jj][pp].pose.view_N_1[1];
		}
		txtOutPut << "\n" << tr("Pose_1_Nz");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "," << KeyPoint_code_serial_update[jj][pp].pose.view_N_1[2];
		}
		txtOutPut << "\n" << tr("Pose_2_Nx");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "," << KeyPoint_code_serial_update[jj][pp].pose.view_N_2[0];
		}
		txtOutPut << "\n" << tr("Pose_2_Ny");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "," << KeyPoint_code_serial_update[jj][pp].pose.view_N_2[1];
		}
		txtOutPut << "\n" << tr("Pose_2_Nz");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "," << KeyPoint_code_serial_update[jj][pp].pose.view_N_2[2];
		}
		txtOutPut << "\n" << tr("Re-Project Error:X");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			if (Result_err_every[jj][pp].x != std::numeric_limits<double>::max())
			{
				txtOutPut << "," << Result_err_every[jj][pp].x;
			}
			else
			{
				txtOutPut << "," << tr("NAN");
			}
		}
		txtOutPut << "\n" << tr("Re-Project Error:Y");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			if (Result_err_every[jj][pp].y != std::numeric_limits<double>::max())
			{
				txtOutPut << "," << Result_err_every[jj][pp].y;
			}
			else
			{
				txtOutPut << "," << tr("NAN");
			}
		}
		txtOutPut << "\n";
		txtOutPut << tr("Rotation(Quaternion)") << "," << Result_R_after_op[jj].x() << ","
			<< Result_R_after_op[jj].y() << "," << Result_R_after_op[jj].z() << "," << Result_R_after_op[jj].w() << "\n";
		txtOutPut << tr("Translation") << "," << Result_T_after_op[jj].x() << ","
			<< Result_T_after_op[jj].y() << "," << Result_T_after_op[jj].z() << "\n";
	}
	txtOutPut << "\n" << "\n" << "\n";


	txtOutPut << tr("Calculate Sequence Index Information") << "\n";
	for (int jj = 0; jj < Keypoints_index_in_queue.size(); jj++)
	{
		txtOutPut << Keypoints_index_in_queue[jj];
		if (jj != (Keypoints_index_in_queue.size() - 1))
		{
			txtOutPut << ",";
		}
	}
	//sbsbsb
	QMessageBox::information(this, tr("Save KeyPoint File"), tr("Save Successfully."), QMessageBox::Ok);
	file.close();
}
void SFM_module::save_SFMResult_file_txt(QString path_save)
{

	QFile file(path_save);
	if (!file.exists())  //文件不存在的时新建
	{
		file.open(QIODevice::WriteOnly | QIODevice::Text);
		QTextStream txtOutPut(&file);
		file.close();
	}
	file.open(QIODevice::WriteOnly | QIODevice::Truncate | QIODevice::Text);
	QTextStream txtOutPut(&file);
	txtOutPut.setRealNumberPrecision(16);
	txtOutPut.setEncoding(QStringConverter::System);

	txtOutPut << tr("Camera Parameter Information") << "\n";
	txtOutPut << tr("IntrinsicMatrix") << "\t" << Result_Camera_K[0] << "\t" << Result_Camera_K[4] << "\t" << Result_Camera_K[2] << "\n";
	txtOutPut << tr("") << "\t" << 0 << "\t" << Result_Camera_K[1] << "\t" << Result_Camera_K[3] << "\n";
	txtOutPut << tr("") << "\t" << 0 << "\t" << 0 << "\t" << 1 << "\n";
	txtOutPut << tr("RadialDistortion") << "\t" << Result_Camera_Dis_K[0] << "\t" << Result_Camera_Dis_K[1] << "\t" << Result_Camera_Dis_K[2] << "\t"
		<< Result_Camera_Dis_K[3] << "\t" << Result_Camera_Dis_K[4] << "\t" << Result_Camera_Dis_K[5]
		<< "\n";
	txtOutPut << tr("TangentialDistortion") << "\t" << Result_Camera_Dis_P[0] << "\t" << Result_Camera_Dis_P[1] << "\n";
	txtOutPut << tr("ThinPrismDistortion") << "\t" << Result_Camera_Dis_T[0] << "\t" << Result_Camera_Dis_T[1]
		<< "\t" << Result_Camera_Dis_T[2] << "\t" << Result_Camera_Dis_T[3] << "\n";

	txtOutPut << "\n" << "\n" << "\n";
	txtOutPut << tr("3D points Information");
	txtOutPut << "\n" << tr("Code");
	for (int pp = 0; pp < Result_xyz_after_op_serial.size(); pp++)
	{
		txtOutPut << "\t" << Result_xyz_after_op_serial[pp].code;
	}
	txtOutPut << "\n" << tr("View Count");
	for (int pp = 0; pp < Result_xyz_after_op_serial.size(); pp++)
	{
		txtOutPut << "\t" << Result_xyz_after_op_serial[pp].weight;
	}
	txtOutPut << "\n" << tr("X");
	for (int pp = 0; pp < Result_xyz_after_op_serial.size(); pp++)
	{
		txtOutPut << "\t" << Result_xyz_after_op_serial[pp].x;
	}
	txtOutPut << "\n" << tr("Y");
	for (int pp = 0; pp < Result_xyz_after_op_serial.size(); pp++)
	{
		txtOutPut << "\t" << Result_xyz_after_op_serial[pp].y;
	}
	txtOutPut << "\n" << tr("Z");
	for (int pp = 0; pp < Result_xyz_after_op_serial.size(); pp++)
	{
		txtOutPut << "\t" << Result_xyz_after_op_serial[pp].z;
	}
	txtOutPut << "\n" << tr("R");
	for (int pp = 0; pp < Result_xyz_after_op_serial.size(); pp++)
	{
		txtOutPut << "\t" << Result_xyz_after_op_serial[pp].R;
	}
	txtOutPut << "\n" << tr("Ori_x");
	for (int pp = 0; pp < Result_xyz_after_op_serial.size(); pp++)
	{
		if (Result_Ori_after_op.size() > pp)
		{
			txtOutPut << "\t" << Result_Ori_after_op[pp].x();
		}
		else
		{
			txtOutPut << "\t" << tr("NAN");
		}
	}
	txtOutPut << "\n" << tr("Ori_y");
	for (int pp = 0; pp < Result_xyz_after_op_serial.size(); pp++)
	{
		if (Result_Ori_after_op.size() > pp)
		{
			txtOutPut << "\t" << Result_Ori_after_op[pp].y();
		}
		else
		{
			txtOutPut << "\t" << tr("NAN");
		}
	}
	txtOutPut << "\n" << tr("Ori_z");
	for (int pp = 0; pp < Result_xyz_after_op_serial.size(); pp++)
	{
		if (Result_Ori_after_op.size() > pp)
		{
			txtOutPut << "\t" << Result_Ori_after_op[pp].z();
		}
		else
		{
			txtOutPut << "\t" << tr("NAN");
		}
	}
	txtOutPut << "\n" << "\n" << "\n";

	txtOutPut << tr("Image-Keypoint Information") << "\n";

	for (int jj = 0; jj < Image_serial_name.size(); jj++)
	{
		txtOutPut << tr("Camera-Index:") << jj + 1 << "\t";
		txtOutPut << Image_serial_name[jj].toLocal8Bit() << "\t";
		txtOutPut << ImageWidth_serial[jj] << "\t";
		txtOutPut << ImageHeight_serial[jj];

		txtOutPut << "\n" << tr("Code");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "\t" << KeyPoint_code_serial_update[jj][pp].code_num;
		}
		txtOutPut << "\n" << tr("X_pixel");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "\t" << KeyPoint_code_serial_update[jj][pp].x;
		}
		txtOutPut << "\n" << tr("Y_pixel");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "\t" << KeyPoint_code_serial_update[jj][pp].y;
		}
		txtOutPut << "\n" << tr("R_a");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "\t" << KeyPoint_code_serial_update[jj][pp].r_a;
		}
		txtOutPut << "\n" << tr("R_b");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "\t" << KeyPoint_code_serial_update[jj][pp].r_b;
		}
		txtOutPut << "\n" << tr("Angle");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "\t" << KeyPoint_code_serial_update[jj][pp].angle_in_pi;
		}
		txtOutPut << "\n" << tr("Pose_1_Nx");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "," << KeyPoint_code_serial_update[jj][pp].pose.view_N_1[0];
		}
		txtOutPut << "\n" << tr("Pose_1_Ny");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "," << KeyPoint_code_serial_update[jj][pp].pose.view_N_1[1];
		}
		txtOutPut << "\n" << tr("Pose_1_Nz");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "," << KeyPoint_code_serial_update[jj][pp].pose.view_N_1[2];
		}
		txtOutPut << "\n" << tr("Pose_2_Nx");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "," << KeyPoint_code_serial_update[jj][pp].pose.view_N_2[0];
		}
		txtOutPut << "\n" << tr("Pose_2_Ny");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "," << KeyPoint_code_serial_update[jj][pp].pose.view_N_2[1];
		}
		txtOutPut << "\n" << tr("Pose_2_Nz");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			txtOutPut << "," << KeyPoint_code_serial_update[jj][pp].pose.view_N_2[2];
		}
		txtOutPut << "\n" << tr("Re-Project Error:X");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			if (Result_err_every[jj][pp].x != std::numeric_limits<double>::max())
			{
				txtOutPut << "\t" << Result_err_every[jj][pp].x;
			}
			else
			{
				txtOutPut << "\t" << tr("NAN");
			}
		}
		txtOutPut << "\n" << tr("Re-Project Error:Y");
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			if (Result_err_every[jj][pp].y != std::numeric_limits<double>::max())
			{
				txtOutPut << "\t" << Result_err_every[jj][pp].y;
			}
			else
			{
				txtOutPut << "\t" << tr("NAN");
			}
		}
		txtOutPut << "\n";
		txtOutPut << tr("Rotation(Quaternion)") << "\t" << Result_R_after_op[jj].x() << "\t"
			<< Result_R_after_op[jj].y() << "\t" << Result_R_after_op[jj].z() << "\t" << Result_R_after_op[jj].w() << "\n";
		txtOutPut << tr("Translation") << "\t" << Result_T_after_op[jj].x() << "\t"
			<< Result_T_after_op[jj].y() << "\t" << Result_T_after_op[jj].z() << "\n";
	}
	txtOutPut << "\n" << "\n" << "\n";


	txtOutPut << tr("Calculate Sequence Index Information") << "\n";
	for (int jj = 0; jj < Keypoints_index_in_queue.size(); jj++)
	{
		txtOutPut << Keypoints_index_in_queue[jj];
		if (jj != (Keypoints_index_in_queue.size() - 1))
		{
			txtOutPut << "\t";
		}
	}
	QMessageBox::information(this, tr("Save Result File"), tr("Save Successfully."), QMessageBox::Ok);
	file.close();
}
void SFM_module::save_SFMResult_file_mat(QString path_save)
{
	MATFile* pmatFile = NULL;
	pmatFile = matOpen(path_save.toLocal8Bit().data(), "w");
	if (pmatFile == NULL)
	{
		QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
		return;
	}
	mxArray* pMxArray_intrinsic = NULL;
	pMxArray_intrinsic = mxCreateDoubleMatrix(3, 3, mxREAL);
	if (pMxArray_intrinsic == NULL)
	{
		QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
		return;
	}
	double* out_data_intrinsic = new double[9];
	out_data_intrinsic[0] = Result_Camera_K[0];
	out_data_intrinsic[3] = Result_Camera_K[4];
	out_data_intrinsic[6] = Result_Camera_K[2];
	out_data_intrinsic[4] = Result_Camera_K[1];
	out_data_intrinsic[7] = Result_Camera_K[3];
	out_data_intrinsic[1] = 0;
	out_data_intrinsic[2] = 0;
	out_data_intrinsic[5] = 0;
	out_data_intrinsic[8] = 1;
	memcpy((void*)(mxGetPr(pMxArray_intrinsic)), (void*)out_data_intrinsic, sizeof(double) * 9);
	memcpy((void*)(mxGetPr(pMxArray_intrinsic)), (void*)out_data_intrinsic, sizeof(double) * 9);
	QString var_mat_name_intrinsic = tr("Camera_IntrinsicMatrix");
	matPutVariable(pmatFile, var_mat_name_intrinsic.toLocal8Bit().toStdString().data(), pMxArray_intrinsic);
	delete[] out_data_intrinsic;
	mxDestroyArray(pMxArray_intrinsic);
	mxArray* pMxArray_radial = NULL;
	pMxArray_radial = mxCreateDoubleMatrix(1, 6, mxREAL);
	if (pMxArray_radial == NULL)
	{
		QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
		return;
	}
	double* out_data_radial = new double[6];
	out_data_radial[0] = Result_Camera_Dis_K[0];
	out_data_radial[1] = Result_Camera_Dis_K[1];
	out_data_radial[2] = Result_Camera_Dis_K[2];
	out_data_radial[3] = Result_Camera_Dis_K[3];
	out_data_radial[4] = Result_Camera_Dis_K[4];
	out_data_radial[5] = Result_Camera_Dis_K[5];
	memcpy((void*)(mxGetPr(pMxArray_radial)), (void*)out_data_radial, sizeof(double) * 6);
	QString var_mat_name_radial = tr("Camera_RadialDistortion");
	matPutVariable(pmatFile, var_mat_name_radial.toLocal8Bit().toStdString().data(), pMxArray_radial);
	delete[] out_data_radial;
	mxDestroyArray(pMxArray_radial);
	mxArray* pMxArray_tangential = NULL;
	pMxArray_tangential = mxCreateDoubleMatrix(1, 2, mxREAL);
	if (pMxArray_tangential == NULL)
	{
		QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
		return;
	}
	double* out_data_tangential = new double[2];
	out_data_tangential[0] = Result_Camera_Dis_P[0];
	out_data_tangential[1] = Result_Camera_Dis_P[1];
	memcpy((void*)(mxGetPr(pMxArray_tangential)), (void*)out_data_tangential, sizeof(double) * 2);
	QString var_mat_name_tangential = tr("Camera_TangentialDistortion");
	matPutVariable(pmatFile, var_mat_name_tangential.toLocal8Bit().toStdString().data(), pMxArray_tangential);
	delete[] out_data_tangential;
	mxDestroyArray(pMxArray_tangential);
	mxArray* pMxArray_thinprism = NULL;
	pMxArray_thinprism = mxCreateDoubleMatrix(1, 4, mxREAL);
	if (pMxArray_thinprism == NULL)
	{
		QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
		return;
	}
	double* out_data_thinprism = new double[4];
	out_data_thinprism[0] = Result_Camera_Dis_T[0];
	out_data_thinprism[1] = Result_Camera_Dis_T[1];
	out_data_thinprism[2] = Result_Camera_Dis_T[2];
	out_data_thinprism[3] = Result_Camera_Dis_T[3];
	memcpy((void*)(mxGetPr(pMxArray_thinprism)), (void*)out_data_thinprism, sizeof(double) * 4);
	QString var_mat_name_thinprism = tr("Camera_ThinPrismDistortion");
	matPutVariable(pmatFile, var_mat_name_thinprism.toLocal8Bit().toStdString().data(), pMxArray_thinprism);
	delete[] out_data_thinprism;
	mxDestroyArray(pMxArray_thinprism);
	for (int jj = 0; jj < Image_serial_name.size(); jj++)
	{
		mxArray* pMxArray = NULL;
		pMxArray = mxCreateDoubleMatrix(8, KeyPoint_code_serial[jj].size(), mxREAL);
		if (pMxArray == NULL)
		{
			QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			return;
		}
		double* out_data = new double[KeyPoint_code_serial[jj].size() * 8];
		for (int pp = 0; pp < KeyPoint_code_serial[jj].size(); pp++)
		{
			out_data[pp * 8] = KeyPoint_code_serial[jj][pp].code_num;
			out_data[pp * 8 + 1] = KeyPoint_code_serial[jj][pp].x;
			out_data[pp * 8 + 2] = KeyPoint_code_serial[jj][pp].y;
			out_data[pp * 8 + 3] = KeyPoint_code_serial[jj][pp].r_a;
			out_data[pp * 8 + 4] = KeyPoint_code_serial[jj][pp].r_b;
			out_data[pp * 8 + 5] = KeyPoint_code_serial[jj][pp].angle_in_pi;
			if (Result_err_every[jj][pp].x == std::numeric_limits<double>::max())
			{
				out_data[pp * 8 + 6] = NAN;
			}
			else
			{
				out_data[pp * 8 + 6] = Result_err_every[jj][pp].x;
			}
			if (Result_err_every[jj][pp].y == std::numeric_limits<double>::max())
			{
				out_data[pp * 8 + 7] = NAN;
			}
			else
			{
				out_data[pp * 8 + 7] = Result_err_every[jj][pp].y;
			}
		}
		memcpy((void*)(mxGetPr(pMxArray)), (void*)out_data, sizeof(double) * KeyPoint_code_serial[jj].size() * 8);
		QString var_mat_name = tr("Image_index") + QString("%1").arg(jj + 1, 8, 10, QLatin1Char('0'));
		matPutVariable(pmatFile, var_mat_name.toLocal8Bit().toStdString().data(), pMxArray);
		delete[] out_data;
		mxDestroyArray(pMxArray);


		mxArray* pMxArray_overall_R = NULL;
		pMxArray_overall_R = mxCreateDoubleMatrix(3, 3, mxREAL);
		if (pMxArray_overall_R == NULL)
		{
			QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			return;
		}
		Eigen::Matrix3d R_odom_curr_tmp = Result_R_after_op[jj].matrix();
		double* out_data_overall_R = new double[9];
		out_data_overall_R[0] = R_odom_curr_tmp(0, 0);
		out_data_overall_R[1] = R_odom_curr_tmp(1, 0);
		out_data_overall_R[2] = R_odom_curr_tmp(2, 0);
		out_data_overall_R[3] = R_odom_curr_tmp(0, 1);
		out_data_overall_R[4] = R_odom_curr_tmp(1, 1);
		out_data_overall_R[5] = R_odom_curr_tmp(2, 1);
		out_data_overall_R[6] = R_odom_curr_tmp(0, 2);
		out_data_overall_R[7] = R_odom_curr_tmp(1, 2);
		out_data_overall_R[8] = R_odom_curr_tmp(2, 2);
		memcpy((void*)(mxGetPr(pMxArray_overall_R)), (void*)out_data_overall_R, sizeof(double) * 9);
		QString var_mat_name_overall_R = tr("Image_index") + 
			QString("%1").arg(jj + 1, 8, 10, QLatin1Char('0')) + tr("_R");
		matPutVariable(pmatFile, var_mat_name_overall_R.toLocal8Bit().toStdString().data(), pMxArray_overall_R);
		delete[] out_data_overall_R;
		mxDestroyArray(pMxArray_overall_R);
		mxArray* pMxArray_overall_T = NULL;
		pMxArray_overall_T = mxCreateDoubleMatrix(1, 3, mxREAL);
		if (pMxArray_overall_T == NULL)
		{
			QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
			return;
		}
		double* out_data_overall_T = new double[3];
		out_data_overall_T[0] = Result_T_after_op[jj].x();
		out_data_overall_T[1] = Result_T_after_op[jj].y();
		out_data_overall_T[2] = Result_T_after_op[jj].z();
		memcpy((void*)(mxGetPr(pMxArray_overall_T)), (void*)out_data_overall_T, sizeof(double) * 3);
		QString var_mat_name_overall_T = tr("Image_index") +
			QString("%1").arg(jj + 1, 8, 10, QLatin1Char('0')) + tr("_T");
		matPutVariable(pmatFile, var_mat_name_overall_T.toLocal8Bit().toStdString().data(), pMxArray_overall_T);
		delete[] out_data_overall_T;
		mxDestroyArray(pMxArray_overall_T);
	}
	mxArray* pMxArray_3d = NULL;
	pMxArray_3d = mxCreateDoubleMatrix(9, Result_xyz_after_op_serial.size(), mxREAL);
	if (pMxArray_3d == NULL)
	{
		QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
		return;
	}
	double* out_data_3d = new double[Result_xyz_after_op_serial.size() * 9];
	for (int pp = 0; pp < Result_xyz_after_op_serial.size(); pp++)
	{
		out_data_3d[pp * 9] = Result_xyz_after_op_serial[pp].code;
		out_data_3d[pp * 9 + 1] = Result_xyz_after_op_serial[pp].weight;
		out_data_3d[pp * 9 + 2] = Result_xyz_after_op_serial[pp].x;
		out_data_3d[pp * 9 + 3] = Result_xyz_after_op_serial[pp].y;
		out_data_3d[pp * 9 + 4] = Result_xyz_after_op_serial[pp].z;
		out_data_3d[pp * 9 + 5] = Result_xyz_after_op_serial[pp].R;

		if (Result_Ori.size() > pp)
		{
			out_data_3d[pp * 9 + 6] = Result_Ori_after_op[pp].x();
			out_data_3d[pp * 9 + 7] = Result_Ori_after_op[pp].y();
			out_data_3d[pp * 9 + 8] = Result_Ori_after_op[pp].z();
		}
		else
		{
			out_data_3d[pp * 9 + 6] = NAN;
			out_data_3d[pp * 9 + 7] = NAN;
			out_data_3d[pp * 9 + 8] = NAN;
		}
	}
	memcpy((void*)(mxGetPr(pMxArray_3d)), (void*)out_data_3d, sizeof(double) * Result_xyz_after_op_serial.size() * 9);
	QString var_mat_name_3d = tr("PointClouds_SFM");
	matPutVariable(pmatFile, var_mat_name_3d.toLocal8Bit().toStdString().data(), pMxArray_3d);
	delete[] out_data_3d;
	mxDestroyArray(pMxArray_3d);
	mxArray* pMxArray_CAL = NULL;
	pMxArray_CAL = mxCreateDoubleMatrix(1, Keypoints_index_in_queue.size(), mxREAL);
	if (pMxArray_CAL == NULL)
	{
		QMessageBox::warning(this, tr("Save Result File"), tr("Save Failed!"), QMessageBox::Ok, QMessageBox::Ok);
		return;
	}
	double* out_data_CAL = new double[Keypoints_index_in_queue.size()];
	for (int pp = 0; pp < Keypoints_index_in_queue.size(); pp++)
	{
		out_data_CAL[pp] = Keypoints_index_in_queue[pp];
	}
	memcpy((void*)(mxGetPr(pMxArray_CAL)), (void*)out_data_CAL, sizeof(double) * Keypoints_index_in_queue.size());
	QString var_mat_name_cal = tr("Calculation_Sequence");
	matPutVariable(pmatFile, var_mat_name_cal.toLocal8Bit().toStdString().data(), pMxArray_CAL);
	delete[] out_data_CAL;
	mxDestroyArray(pMxArray_CAL);
	QMessageBox::information(this, tr("Save Result File"), tr("Save Successfully."), QMessageBox::Ok);
	matClose(pmatFile);
}

void SFM_module::export_inf_file()
{
	if (ui.save_keypoint_checkBox->isChecked())
	{
		select_path_save_keypoint();
	}
	if (ui.save_calculation_result_checkBox->isChecked())
	{
		select_path_save_result();
	}
}

void SFM_module::select_path_save_keypoint()
{
	switch (ui.save_type_comboBox->currentIndex())
	{
	case 0:
	{
		QString file_for_msdata = QFileDialog::getSaveFileName(this, tr("Save KeyPoints data"), tr(""), "MSDATA(*.msdata)");
		if (file_for_msdata.isEmpty())
		{
			return;
		}
		save_KeyPoint_file(file_for_msdata);
	}
	break;
	case 1:
	{
		QString file_for_msdata = QFileDialog::getSaveFileName(this, tr("Save KeyPoints data"), tr(""), "CSV(*.csv)");
		if (file_for_msdata.isEmpty())
		{
			return;
		}
		save_KeyPoint_file_csv(file_for_msdata);
	}
	break;
	case 2:
	{
		QString file_for_msdata = QFileDialog::getSaveFileName(this, tr("Save KeyPoints data"), tr(""), "TXT(*.txt)");
		if (file_for_msdata.isEmpty())
		{
			return;
		}
		save_KeyPoint_file_txt(file_for_msdata);
	}
	break;
	case 3:
	{
		QString file_for_msdata = QFileDialog::getSaveFileName(this, tr("Save KeyPoints data"), tr(""), "MAT(*.mat)");
		if (file_for_msdata.isEmpty())
		{
			return;
		}
		save_KeyPoint_file_mat(file_for_msdata);
	}
	break;
	default:
	{
		QString file_for_msdata = QFileDialog::getSaveFileName(this, tr("Save KeyPoints data"), tr(""), "MSDATA(*.msdata)");
		if (file_for_msdata.isEmpty())
		{
			return;
		}
		save_KeyPoint_file(file_for_msdata);
	}
	break;
	}
}
void SFM_module::select_path_read_keypoint()
{
	QString file_for_msdata = QFileDialog::getOpenFileName(this, tr("Select KeyPoints data"), tr(""), "MSDATA(*.msdata)");
	if (file_for_msdata.isEmpty())
	{
		return;
	}
	clear_all_data();
	read_KeyPoint_file(file_for_msdata);
}
void SFM_module::select_path_save_result()
{
	if (Keypoints_index_in_queue.size() == 0)
	{
		QMessageBox::warning(this, tr("Save Result"), tr("No calculation sequence, result file will not be generated!"), QMessageBox::Ok, QMessageBox::Ok);
		return;
	}
	switch (ui.save_type_comboBox->currentIndex())
	{
	case 0:
	{
		QString file_for_msdata = QFileDialog::getSaveFileName(this, tr("Save Result data"), tr(""), "MSDATA(*.msdata)");
		if (file_for_msdata.isEmpty())
		{
			return;
		}
		save_SFMResult_file(file_for_msdata);
	}
	break;
	case 1:
	{
		QString file_for_msdata = QFileDialog::getSaveFileName(this, tr("Save Result data"), tr(""), "CSV(*.csv)");
		if (file_for_msdata.isEmpty())
		{
			return;
		}
		save_SFMResult_file_csv(file_for_msdata);
	}
	break;
	case 2:
	{
		QString file_for_msdata = QFileDialog::getSaveFileName(this, tr("Save Result data"), tr(""), "TXT(*.txt)");
		if (file_for_msdata.isEmpty())
		{
			return;
		}
		save_SFMResult_file_txt(file_for_msdata);
	}
	break;
	case 3:
	{
		QString file_for_msdata = QFileDialog::getSaveFileName(this, tr("Save Result data"), tr(""), "MAT(*.mat)");
		if (file_for_msdata.isEmpty())
		{
			return;
		}
		save_SFMResult_file_mat(file_for_msdata);
	}
	break;
	default:
	{
		QString file_for_msdata = QFileDialog::getSaveFileName(this, tr("Save Result data"), tr(""), "MSDATA(*.msdata)");
		if (file_for_msdata.isEmpty())
		{
			return;
		}
		save_SFMResult_file(file_for_msdata);
	}
	break;
	}
}
void SFM_module::select_path_read_result()
{
	QString file_for_msdata = QFileDialog::getOpenFileName(this, tr("Select KeyPoints data"), tr(""), "MSDATA(*.msdata)");
	if (file_for_msdata.isEmpty())
	{
		return;
	}
	clear_all_data();
	read_SFMResult_file(file_for_msdata);
}

deltect_coded_circle_thread::deltect_coded_circle_thread(QString image_path,
	float ratio_k, float ratio_k1, float ratio_k2,
	float min_radius, float max_radius, float ellipse_error_pixel,
	MarkPointColorType color_type, CodePointBitesType code_bites_type,
	DetectContoursMethod image_process_method, SubPixelPosMethod subpixel_pos_method
	, double max_aspect_ratio, int min_points, int min_contour_num,
	float delta_Mt, float fore_stdDev, float back_stdDev)
{
	image_path_thread = image_path;
	ratio_k_thread = ratio_k;
	ratio_k1_thread = ratio_k1;
	ratio_k2_thread = ratio_k2;
	min_radius_thread = min_radius;
	max_radius_thread = max_radius;
	ellipse_error_pixel_thread = ellipse_error_pixel;
	color_type_thread = color_type;
	code_bites_type_thread = code_bites_type;
	image_process_method_thread = image_process_method;
	subpixel_pos_method_thread = subpixel_pos_method;
	max_aspect_ratio_thread = max_aspect_ratio;
	min_points_thread = min_points;
	min_contour_num_thread = min_contour_num;
	delta_Mt_thread = delta_Mt;
	fore_stdDev_thread = fore_stdDev;
	back_stdDev_thread = back_stdDev;
}

void deltect_coded_circle_thread::run()
{
	ori_image_thread = cv::imread(image_path_thread.toLocal8Bit().toStdString());
	success_cal = ImageDetectMethod::detectcodecircle(ori_image_thread, code_point_mat, contours_pixel,
		ratio_k_thread, ratio_k1_thread, ratio_k2_thread,
		min_radius_thread, max_radius_thread,
		ellipse_error_pixel_thread, color_type_thread,
		code_bites_type_thread, image_process_method_thread, subpixel_pos_method_thread, max_aspect_ratio_thread,
		min_points_thread, min_contour_num_thread, delta_Mt_thread, fore_stdDev_thread, back_stdDev_thread);
	is_finish = true;
}

optimize_SFM_incremental_thread::optimize_SFM_incremental_thread(std::vector<std::vector<Coded_detect_inf>> code_circle_serial,
	std::vector<sfm_3d_group> point_clouds,
	std::vector<Eigen::Quaterniond> camQvec, std::vector<Eigen::Vector3d> camTvec
	, std::vector<std::vector<cv::Point2d>> Re_project_Map
	, double* camK, double* camDK, double* camDP, double* camDT
	, ceres::Solver::Summary* summary
	, unsigned int max_iter_num
	, double stop_value
	, unsigned int num_threads
	, unsigned int timeout
	, unsigned char Loss_type
	, double Loss_value
	, unsigned int Dis_K_num
	, unsigned int Dis_P_num
	, unsigned int Dis_T_num
	, bool Fixed_all
	, bool Use_F
	, bool Use_Cx_Cy
	, bool Use_shear
	, bool Use_same_F
	, double* F_range
	, double* C_range
	, bool is_part)
{
	code_circle_serial_thread = code_circle_serial;
	point_clouds_thread = point_clouds;
	camQvec_thread = camQvec;
	camTvec_thread = camTvec;
	Re_project_Map_thread = Re_project_Map;
	camK_thread = camK;
	camDK_thread = camDK;
	camDP_thread = camDP;
	camDT_thread = camDT;
	summary_thread = summary;
	max_iter_num_thread = max_iter_num;
	stop_value_thread = stop_value;
	num_threads_thread = num_threads;
	timeout_thread = timeout;
	Loss_type_thread = Loss_type;
	Loss_value_thread = Loss_value;
	Dis_K_num_thread = Dis_K_num;
	Dis_P_num_thread = Dis_P_num;
	Dis_T_num_thread = Dis_T_num;
	Use_F_thread = Use_F;
	Use_Cx_Cy_thread = Use_Cx_Cy;
	Use_shear_thread = Use_shear;
	Use_same_F_thread = Use_same_F;
	F_range_thread = F_range;
	C_range_thread = C_range;
	Fixed_all_trhead = Fixed_all;
	is_part_thread = is_part;
}


void optimize_SFM_incremental_thread::run()
{
	//SBSBSB
	if (is_part_thread)
	{
		SFM_optimize::calculate_SFM_part(code_circle_serial_thread,
			point_clouds_thread,
			camQvec_thread, camTvec_thread
			, Re_project_Map_thread
			, camK_thread, camDK_thread, camDP_thread, camDT_thread
			, summary_thread
			, max_iter_num_thread
			, stop_value_thread
			, num_threads_thread
			, timeout_thread
			, Loss_type_thread
			, Loss_value_thread
			, Dis_K_num_thread
			, Dis_P_num_thread
			, Dis_T_num_thread
			, Fixed_all_trhead
			, Use_F_thread
			, Use_Cx_Cy_thread
			, Use_shear_thread
			, Use_same_F_thread
			, F_range_thread
			, C_range_thread);
	}
	else
	{
		SFM_optimize::calculate_SFM(code_circle_serial_thread,
			point_clouds_thread,
			camQvec_thread, camTvec_thread
			, Re_project_Map_thread
			, camK_thread, camDK_thread, camDP_thread, camDT_thread
			, summary_thread
			, max_iter_num_thread
			, stop_value_thread
			, num_threads_thread
			, timeout_thread
			, Loss_type_thread
			, Loss_value_thread
			, Dis_K_num_thread
			, Dis_P_num_thread
			, Dis_T_num_thread
			, Fixed_all_trhead
			, Use_F_thread
			, Use_Cx_Cy_thread
			, Use_shear_thread
			, Use_same_F_thread
			, F_range_thread
			, C_range_thread);
	}
	is_finish = true;
}

calculate_circle_thread::calculate_circle_thread(std::vector<sfm_3d_group> point_cloudss, std::vector<std::vector<std::vector<cv::Point>>> contour_circle_serials,
	std::vector<std::vector<Coded_detect_inf>> code_circle_serials, std::vector<Eigen::Vector3d> ori_serials,
	std::vector<Eigen::Quaterniond> camQvecs, std::vector<Eigen::Vector3d> camTvecs
	, double* camKs, double* camDKs, double* camDPs, double* camDTs
	, ceres::Solver::Summary* summarys
	, unsigned int max_iter_nums
	, double stop_values
	, unsigned int num_threadss
	, unsigned int timeouts
	, unsigned char Loss_types
	, double Loss_values)
{
	point_clouds = point_cloudss;
	contour_circle_serial = contour_circle_serials;
	code_circle_serial = code_circle_serials;
	ori_serial = ori_serials;
	camQvec = camQvecs;
	camTvec = camTvecs;
	camK = camKs;
	camDK = camDKs;
	camDP = camDPs;
	camDT = camDTs;
	summary = summarys;
	max_iter_num = max_iter_nums;
	stop_value = stop_values;
	num_threads = num_threadss;
	timeout = timeouts;
	Loss_type = Loss_types;
	Loss_value = Loss_values;
}

void calculate_circle_thread::run()
{
	SFM_optimize::calculate_Circle(point_clouds, contour_circle_serial,
		code_circle_serial, ori_serial,
		camQvec, camTvec
		, camK, camDK, camDP, camDT
		, summary
		, max_iter_num
		, stop_value
		, num_threads
		, timeout
		, Loss_type
		, Loss_value);
	is_finish = true;
}

update_coded_circle_pixel_thread::update_coded_circle_pixel_thread(Eigen::Matrix<double, 3, 4> Camera_K
	, Eigen::Matrix<double, 3, 3> Camera_R
	, Eigen::Vector3d Camera_T, std::vector<Coded_detect_inf> key_points
	, std::vector<Eigen::Vector3d> ori_vec, std::vector<sfm_3d_group> point_3d, std::vector<std::vector<cv::Point>> contour,
	QString image_path,
	float ratio_k, float ratio_k1, float ratio_k2,
	float min_radius, float max_radius, float ellipse_error_pixel,
	MarkPointColorType color_type, CodePointBitesType code_bites_type,
	DetectContoursMethod image_process_method, SubPixelPosMethod subpixel_pos_method
	, double max_aspect_ratio, int min_points, int min_contour_num,
	float delta_Mt, float fore_stdDev, float back_stdDev)
{
	Camera_K_thread = Camera_K;
	Camera_R_thread = Camera_R;
	Camera_T_thread = Camera_T;
	key_points_thread = key_points;
	ori_vec_thread = ori_vec;
	point_3d_thread = point_3d;
	contour_thread = contour;
	image_path_thread = image_path;
	ratio_k_thread = ratio_k;
	ratio_k1_thread = ratio_k1;
	ratio_k2_thread = ratio_k2;
	min_radius_thread = min_radius;
	max_radius_thread = max_radius;
	ellipse_error_pixel_thread = ellipse_error_pixel;
	color_type_thread = color_type;
	code_bites_type_thread = code_bites_type;
	image_process_method_thread = image_process_method;
	subpixel_pos_method_thread = subpixel_pos_method;
	max_aspect_ratio_thread = max_aspect_ratio;
	min_points_thread = min_points;
	min_contour_num_thread = min_contour_num;
	delta_Mt_thread = delta_Mt;
	fore_stdDev_thread = fore_stdDev;
	back_stdDev_thread = back_stdDev;
}
void update_coded_circle_pixel_thread::run()
{
	QString img_name = image_path_thread;
	cv::Mat image_now = cv::imread(img_name.toLocal8Bit().toStdString());
	if (!image_now.data)
	{
		is_finish = true;
		return;
	}
	if (image_now.channels() == 3)
	{
		cv::cvtColor(image_now, image_now, CV_BGR2GRAY);
	}
	for (int jj = 0; jj < key_points_thread.size(); jj++)
	{
		int find_sfm_flag = -1;
		for (int ss = 0; ss < point_3d_thread.size(); ss++)
		{
			if (key_points_thread[jj].code_num == point_3d_thread[ss].code)
			{
				find_sfm_flag = ss;
				break;
			}
		}
		if (find_sfm_flag == -1)
		{
			continue;
		}

		double min_x = std::numeric_limits<double>::max();
		double min_y = std::numeric_limits<double>::max();
		double max_x = -std::numeric_limits<double>::max();
		double max_y = -std::numeric_limits<double>::max();
		for (int pp = 0; pp < contour_thread[jj].size(); pp++)
		{
			min_x = min_x < contour_thread[jj][pp].x ? min_x : contour_thread[jj][pp].x;
			min_y = min_y < contour_thread[jj][pp].y ? min_y : contour_thread[jj][pp].y;
			max_x = max_x > contour_thread[jj][pp].x ? max_x : contour_thread[jj][pp].x;
			max_y = max_y > contour_thread[jj][pp].y ? max_y : contour_thread[jj][pp].y;
		}
		double muli_cof = 1.25;
		min_x = key_points_thread[jj].x - muli_cof * ratio_k2_thread * (key_points_thread[jj].x - min_x);
		max_x = key_points_thread[jj].x + muli_cof * ratio_k2_thread * (-key_points_thread[jj].x + max_x);
		min_y = key_points_thread[jj].y - muli_cof * ratio_k2_thread * (key_points_thread[jj].y - min_y);
		max_y = key_points_thread[jj].y + muli_cof * ratio_k2_thread * (-key_points_thread[jj].y + max_y);
		min_x = floor(min_x);
		min_y = floor(min_y);
		max_x = ceil(max_x);
		max_y = ceil(max_y);
		min_x = min_x < 0 ? 0 : min_x;
		min_y = min_y < 0 ? 0 : min_y;
		max_x = max_x >= image_now.cols ? (image_now.cols - 1) : max_x;
		max_y = max_y >= image_now.rows ? (image_now.rows - 1) : max_y;
		cv::Mat Circle_image_Old = image_now(cv::Range(min_y, max_y), cv::Range(min_x, max_x));

		int window_R_new = ceil(muli_cof * ((max_y - min_y) > (max_x - min_x) ? (max_y - min_y) : (max_x - min_x)) / 2.0);
		//Eigen::MatrixXd new_ori_temp = camQvec[ii].toRotationMatrix() * Result_Ori[ii].matrix();
		Eigen::MatrixXd new_ori_temp = ori_vec_thread[find_sfm_flag].matrix();
		Eigen::Vector3d new_ori(new_ori_temp(0, 0), new_ori_temp(1, 0), new_ori_temp(2, 0));
		Eigen::Vector3d cros_i(1, 0, 0);
		Eigen::Vector3d cros_j(0, 1, 0);
		Eigen::Vector3d cros_a, cros_b;
		cros_a = new_ori.cross(cros_i);
		if (cros_a.x() == 0 && cros_a.y() == 0 && cros_a.z() == 0)
		{
			cros_a = new_ori.cross(cros_j);
		}
		cros_b = new_ori.cross(cros_a);
		cros_a.normalize();
		cros_b.normalize();

		Eigen::MatrixXd new_sfm_xyz = Eigen::Vector3d(point_3d_thread[find_sfm_flag].x, point_3d_thread[find_sfm_flag].y, point_3d_thread[find_sfm_flag].z).matrix();
		Eigen::Vector4d old_point_1_4d(new_sfm_xyz(0, 0) + muli_cof * ratio_k2_thread * point_3d_thread[find_sfm_flag].R * cros_a(0),
			new_sfm_xyz(1, 0) + muli_cof * ratio_k2_thread * point_3d_thread[find_sfm_flag].R * cros_a(1), new_sfm_xyz(2, 0) + muli_cof * ratio_k2_thread * point_3d_thread[find_sfm_flag].R * cros_a(2), 1);
		Eigen::Vector4d old_point_2_4d(new_sfm_xyz(0, 0) + muli_cof * ratio_k2_thread * point_3d_thread[find_sfm_flag].R * cros_b(0),
			new_sfm_xyz(1, 0) + muli_cof * ratio_k2_thread * point_3d_thread[find_sfm_flag].R * cros_b(1), new_sfm_xyz(2, 0) + muli_cof * ratio_k2_thread * point_3d_thread[find_sfm_flag].R * cros_b(2), 1);
		Eigen::Vector4d old_point_3_4d(new_sfm_xyz(0, 0) - muli_cof * ratio_k2_thread * point_3d_thread[find_sfm_flag].R * cros_a(0),
			new_sfm_xyz(1, 0) - muli_cof * ratio_k2_thread * point_3d_thread[find_sfm_flag].R * cros_a(1), new_sfm_xyz(2, 0) - muli_cof * ratio_k2_thread * point_3d_thread[find_sfm_flag].R * cros_a(2), 1);
		Eigen::Vector4d old_point_4_4d(new_sfm_xyz(0, 0) - muli_cof * ratio_k2_thread * point_3d_thread[find_sfm_flag].R * cros_b(0),
			new_sfm_xyz(1, 0) - muli_cof * ratio_k2_thread * point_3d_thread[find_sfm_flag].R * cros_b(1), new_sfm_xyz(2, 0) - muli_cof * ratio_k2_thread * point_3d_thread[find_sfm_flag].R * cros_b(2), 1);

		Eigen::Vector3d old_point_1_3d = Camera_K_thread * old_point_1_4d.matrix();
		Eigen::Vector3d old_point_2_3d = Camera_K_thread * old_point_2_4d.matrix();
		Eigen::Vector3d old_point_3_3d = Camera_K_thread * old_point_3_4d.matrix();
		Eigen::Vector3d old_point_4_3d = Camera_K_thread * old_point_4_4d.matrix();
		cv::Point2f old_point_1(old_point_1_3d.x() / old_point_1_3d.z() - min_x, old_point_1_3d.y() / old_point_1_3d.z() - min_y);
		cv::Point2f old_point_2(old_point_2_3d.x() / old_point_2_3d.z() - min_x, old_point_2_3d.y() / old_point_2_3d.z() - min_y);
		cv::Point2f old_point_3(old_point_3_3d.x() / old_point_3_3d.z() - min_x, old_point_3_3d.y() / old_point_3_3d.z() - min_y);
		cv::Point2f old_point_4(old_point_4_3d.x() / old_point_4_3d.z() - min_x, old_point_4_3d.y() / old_point_4_3d.z() - min_y);
		cv::Point2f new_point_1(window_R_new, 0);
		cv::Point2f new_point_2(0, window_R_new);
		cv::Point2f new_point_3(window_R_new, 2 * window_R_new);
		cv::Point2f new_point_4(2 * window_R_new, window_R_new);

		cv::Point2f AffinePointsSrc[4] = { old_point_1, old_point_2, old_point_3, old_point_4 };
		cv::Point2f AffinePointsDst[4] = { new_point_1, new_point_2, new_point_3, new_point_4 };
		cv::Mat Circle_image_New;
		cv::Mat per_tf = cv::getPerspectiveTransform(AffinePointsSrc, AffinePointsDst);
		cv::warpPerspective(Circle_image_Old, Circle_image_New, per_tf, cv::Size(2 * window_R_new + 1, 2 * window_R_new + 1), CV_INTER_AREA);

		cv::Mat code_point_temp;
	    std::vector<std::vector<cv::Point>> pixels;
		ImageDetectMethod::detectcodecircle(Circle_image_New, code_point_temp, pixels, ratio_k_thread,
			ratio_k1_thread, ratio_k2_thread,
			2, window_R_new * 2 + 1,
			ellipse_error_pixel_thread, color_type_thread, code_bites_type_thread,
			image_process_method_thread, subpixel_pos_method_thread, max_aspect_ratio_thread
			, min_points_thread, min_contour_num_thread
			, delta_Mt_thread
			, fore_stdDev_thread, back_stdDev_thread);
		bool can_be_ex = false;
		for (unsigned int gg = 0; gg < code_point_temp.rows; gg++)
		{
			//if (code_point_temp.at<float>(gg, 0) == 9)
			//{
			//	cv::namedWindow("Old", cv::WINDOW_NORMAL);
			//	cv::imshow("Old", Circle_image_Old);
			//	cv::namedWindow("New", cv::WINDOW_NORMAL);
			//	cv::imshow("New", Circle_image_New);
			//	cv::waitKey(0);
			//}
			if (code_point_temp.at<float>(gg, 0) == key_points_thread[jj].code_num)
			{
				std::cout << code_point_temp.at<float>(gg, 0) << "\t" << key_points_thread[jj].r_a << "\t" << key_points_thread[jj].r_b << "\t" << code_point_temp.at<float>(gg, 4) << "\t" << code_point_temp.at<float>(gg, 5) << std::endl;
				can_be_ex = true;
				cv::Mat_<double> mat_pt(3, 1);
				mat_pt(0, 0) = code_point_temp.at<float>(gg, 1);
				mat_pt(1, 0) = code_point_temp.at<float>(gg, 2);
				mat_pt(2, 0) = 1;
				cv::Mat mat_pt_view = per_tf.inv() * mat_pt;
				double a1 = mat_pt_view.at<double>(0, 0);
				double a2 = mat_pt_view.at<double>(1, 0);
				double a3 = mat_pt_view.at<double>(2, 0);
				double new_x = a1 * 1.0 / a3 + min_x;
				double new_y = a2 * 1.0 / a3 + min_y;
				key_points_thread[jj].x = new_x;
				key_points_thread[jj].y = new_y;
				break;
			}
		}
	}
	is_finish = true;
}