#include "SFM_optimize.h"
#include <core/eigen.hpp>
void SFM_optimize::calculate_Circle(std::vector<sfm_3d_group> &point_clouds, std::vector<std::vector<std::vector<cv::Point>>> contour_circle_serial,
	std::vector<std::vector<Coded_detect_inf>> code_circle_serial, std::vector<Eigen::Vector3d>& ori_serial,
	std::vector<Eigen::Quaterniond> camQvec, std::vector<Eigen::Vector3d> camTvec
	, double* camK, double* camDK, double* camDP, double* camDT
	, ceres::Solver::Summary* summary
	, unsigned int max_iter_num
	, double stop_value
	, unsigned int num_threads
	, unsigned int timeout
	, unsigned char Loss_type
	, double Loss_value)
{
	double* Circle_params = new double[point_clouds.size() * 6];
	ceres::Problem problem;
	std::vector<Eigen::Matrix<double, 3, 4>> K_serial;
	Eigen::Matrix<double, 3, 3> I_matrix = Eigen::Matrix<double, 3, 3>::Zero();
	I_matrix(0, 0) = camK[0];
	I_matrix(0, 1) = camK[4];
	I_matrix(0, 2) = camK[2];
	I_matrix(1, 1) = camK[1];
	I_matrix(1, 2) = camK[3];
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
	for (int ii = 0; ii < point_clouds.size(); ii++)
	{
		double mod = sqrt(pow(ori_serial[ii].x(), 2) + pow(ori_serial[ii].y(), 2) + pow(ori_serial[ii].z(), 2));
		ori_serial[ii].x() /= mod;
		ori_serial[ii].y() /= mod;
		ori_serial[ii].z() /= mod;
	}
	std::vector<std::vector<std::vector<cv::Point2d>>> contour_circle_double_serial;
	for (int pp = 0; pp < contour_circle_serial.size(); pp++)
	{
		std::vector<std::vector<cv::Point2d>> double_serial_view = {};
		for (int qq = 0; qq < contour_circle_serial[pp].size(); qq++)
		{
			std::vector<cv::Point2d> double_serial = {};
			for (int tt = 0; tt < contour_circle_serial[pp][qq].size(); tt++)
			{
				double_serial.push_back(cv::Point2d(contour_circle_serial[pp][qq][tt].x, contour_circle_serial[pp][qq][tt].y));
			}
			double_serial_view.push_back(double_serial);
		}
		contour_circle_double_serial.push_back(double_serial_view);
	}
	cv::Mat K_camera = cv::Mat(cv::Matx33d(
		camK[0], camK[4], camK[2],
		0, camK[1], camK[3],
		0, 0, 1));
	cv::Mat dis_coff = cv::Mat::eye(1, 12, CV_64FC1);
	dis_coff.at<double>(0) = camDK[0];
	dis_coff.at<double>(1) = camDK[1];
	dis_coff.at<double>(2) = camDP[0];
	dis_coff.at<double>(3) = camDP[1];
	dis_coff.at<double>(4) = camDK[2];
	dis_coff.at<double>(5) = camDK[3];
	dis_coff.at<double>(6) = camDK[4];
	dis_coff.at<double>(7) = camDK[5];
	dis_coff.at<double>(8) = camDT[0];
	dis_coff.at<double>(9) = camDT[2];
	dis_coff.at<double>(10) = camDT[1];
	dis_coff.at<double>(11) = camDT[3];

	for (int pp = 0; pp < contour_circle_double_serial.size(); pp++)
	{
		for (int qq = 0; qq < contour_circle_double_serial[pp].size(); qq++)
		{
			cv::undistortPoints(contour_circle_double_serial[pp][qq], contour_circle_double_serial[pp][qq], K_camera, dis_coff, K_camera);
		}
	}

	for (int ii = 0; ii < point_clouds.size(); ii++)
	{
		Circle_params[ii * 6 + 0] = point_clouds[ii].x;
		Circle_params[ii * 6 + 1] = point_clouds[ii].y;
		Circle_params[ii * 6 + 2] = point_clouds[ii].z;
		Circle_params[ii * 6 + 4] = asin(ori_serial[ii].x());
		Circle_params[ii * 6 + 5] = atan2(ori_serial[ii].y(), ori_serial[ii].z());
	}
	std::vector<double> circle_R;
	std::vector<double> circle_R_Number;
	for (int ii = 0; ii < point_clouds.size(); ii++)
	{
		circle_R.push_back(0);
		circle_R_Number.push_back(0);
	}
	for (int ii = 0; ii < point_clouds.size(); ii++)
	{
		for (int pp = 0; pp < code_circle_serial.size(); pp++)
		{
			for (int qq = 0; qq < code_circle_serial[pp].size(); qq++)
			{
				if (code_circle_serial[pp][qq].code_num == point_clouds[ii].code)
				{
					for (int tt = 0; tt < contour_circle_double_serial[pp][qq].size(); tt++)
					{
						double N1 = sin(Circle_params[ii * 6 + 4]);
						double N2 = cos(Circle_params[ii * 6 + 4]) * sin(Circle_params[ii * 6 + 5]);
						double N3 = cos(Circle_params[ii * 6 + 4]) * cos(Circle_params[ii * 6 + 5]);
						double X0 = contour_circle_double_serial[pp][qq][tt].x;
						double Y0 = contour_circle_double_serial[pp][qq][tt].y;
						double x_w = (K_serial[pp](0, 1) * K_serial[pp](1, 3) * N3 - K_serial[pp](0, 2) * K_serial[pp](1, 3) * N2 - K_serial[pp](0, 3) * K_serial[pp](1, 1) * N3 + K_serial[pp](0, 3) * K_serial[pp](1, 2) * N2 + K_serial[pp](1, 1) * K_serial[pp](2, 3) * N3 * X0 - K_serial[pp](1, 2) * K_serial[pp](2, 3) * N2 * X0 - K_serial[pp](1, 3) * K_serial[pp](2, 1) * N3 * X0 + K_serial[pp](1, 3) * K_serial[pp](2, 2) * N2 * X0 + K_serial[pp](0, 1) * K_serial[pp](1, 2) * N1 * Circle_params[ii * 6 + 0] - K_serial[pp](0, 2) * K_serial[pp](1, 1) * N1 * Circle_params[ii * 6 + 0] - K_serial[pp](0, 1) * K_serial[pp](2, 3) * N3 * Y0 + K_serial[pp](0, 2) * K_serial[pp](2, 3) * N2 * Y0 + K_serial[pp](0, 3) * K_serial[pp](2, 1) * N3 * Y0 - K_serial[pp](0, 3) * K_serial[pp](2, 2) * N2 * Y0 + K_serial[pp](0, 1) * K_serial[pp](1, 2) * N2 * Circle_params[ii * 6 + 1] - K_serial[pp](0, 2) * K_serial[pp](1, 1) * N2 * Circle_params[ii * 6 + 1] + K_serial[pp](0, 1) * K_serial[pp](1, 2) * N3 * Circle_params[ii * 6 + 2] - K_serial[pp](0, 2) * K_serial[pp](1, 1) * N3 * Circle_params[ii * 6 + 2] + K_serial[pp](1, 1) * K_serial[pp](2, 2) * N1 * X0 * Circle_params[ii * 6 + 0] - K_serial[pp](1, 2) * K_serial[pp](2, 1) * N1 * X0 * Circle_params[ii * 6 + 0] - K_serial[pp](0, 1) * K_serial[pp](2, 2) * N1 * Circle_params[ii * 6 + 0] * Y0 + K_serial[pp](0, 2) * K_serial[pp](2, 1) * N1 * Circle_params[ii * 6 + 0] * Y0 + K_serial[pp](1, 1) * K_serial[pp](2, 2) * N2 * X0 * Circle_params[ii * 6 + 1] - K_serial[pp](1, 2) * K_serial[pp](2, 1) * N2 * X0 * Circle_params[ii * 6 + 1] - K_serial[pp](0, 1) * K_serial[pp](2, 2) * N2 * Y0 * Circle_params[ii * 6 + 1] + K_serial[pp](0, 2) * K_serial[pp](2, 1) * N2 * Y0 * Circle_params[ii * 6 + 1] + K_serial[pp](1, 1) * K_serial[pp](2, 2) * N3 * X0 * Circle_params[ii * 6 + 2] - K_serial[pp](1, 2) * K_serial[pp](2, 1) * N3 * X0 * Circle_params[ii * 6 + 2] - K_serial[pp](0, 1) * K_serial[pp](2, 2) * N3 * Y0 * Circle_params[ii * 6 + 2] + K_serial[pp](0, 2) * K_serial[pp](2, 1) * N3 * Y0 * Circle_params[ii * 6 + 2]) / (K_serial[pp](0, 0) * K_serial[pp](1, 1) * N3 - K_serial[pp](0, 0) * K_serial[pp](1, 2) * N2 - K_serial[pp](0, 1) * K_serial[pp](1, 0) * N3 + K_serial[pp](0, 1) * K_serial[pp](1, 2) * N1 + K_serial[pp](0, 2) * K_serial[pp](1, 0) * N2 - K_serial[pp](0, 2) * K_serial[pp](1, 1) * N1 + K_serial[pp](1, 0) * K_serial[pp](2, 1) * N3 * X0 - K_serial[pp](1, 0) * K_serial[pp](2, 2) * N2 * X0 - K_serial[pp](1, 1) * K_serial[pp](2, 0) * N3 * X0 + K_serial[pp](1, 1) * K_serial[pp](2, 2) * N1 * X0 + K_serial[pp](1, 2) * K_serial[pp](2, 0) * N2 * X0 - K_serial[pp](1, 2) * K_serial[pp](2, 1) * N1 * X0 - K_serial[pp](0, 0) * K_serial[pp](2, 1) * N3 * Y0 + K_serial[pp](0, 0) * K_serial[pp](2, 2) * N2 * Y0 + K_serial[pp](0, 1) * K_serial[pp](2, 0) * N3 * Y0 - K_serial[pp](0, 1) * K_serial[pp](2, 2) * N1 * Y0 - K_serial[pp](0, 2) * K_serial[pp](2, 0) * N2 * Y0 + K_serial[pp](0, 2) * K_serial[pp](2, 1) * N1 * Y0);
						double y_w = -(K_serial[pp](0, 0) * K_serial[pp](1, 3) * N3 - K_serial[pp](0, 2) * K_serial[pp](1, 3) * N1 - K_serial[pp](0, 3) * K_serial[pp](1, 0) * N3 + K_serial[pp](0, 3) * K_serial[pp](1, 2) * N1 + K_serial[pp](1, 0) * K_serial[pp](2, 3) * N3 * X0 - K_serial[pp](1, 2) * K_serial[pp](2, 3) * N1 * X0 - K_serial[pp](1, 3) * K_serial[pp](2, 0) * N3 * X0 + K_serial[pp](1, 3) * K_serial[pp](2, 2) * N1 * X0 + K_serial[pp](0, 0) * K_serial[pp](1, 2) * N1 * Circle_params[ii * 6 + 0] - K_serial[pp](0, 2) * K_serial[pp](1, 0) * N1 * Circle_params[ii * 6 + 0] - K_serial[pp](0, 0) * K_serial[pp](2, 3) * N3 * Y0 + K_serial[pp](0, 2) * K_serial[pp](2, 3) * N1 * Y0 + K_serial[pp](0, 3) * K_serial[pp](2, 0) * N3 * Y0 - K_serial[pp](0, 3) * K_serial[pp](2, 2) * N1 * Y0 + K_serial[pp](0, 0) * K_serial[pp](1, 2) * N2 * Circle_params[ii * 6 + 1] - K_serial[pp](0, 2) * K_serial[pp](1, 0) * N2 * Circle_params[ii * 6 + 1] + K_serial[pp](0, 0) * K_serial[pp](1, 2) * N3 * Circle_params[ii * 6 + 2] - K_serial[pp](0, 2) * K_serial[pp](1, 0) * N3 * Circle_params[ii * 6 + 2] + K_serial[pp](1, 0) * K_serial[pp](2, 2) * N1 * X0 * Circle_params[ii * 6 + 0] - K_serial[pp](1, 2) * K_serial[pp](2, 0) * N1 * X0 * Circle_params[ii * 6 + 0] - K_serial[pp](0, 0) * K_serial[pp](2, 2) * N1 * Circle_params[ii * 6 + 0] * Y0 + K_serial[pp](0, 2) * K_serial[pp](2, 0) * N1 * Circle_params[ii * 6 + 0] * Y0 + K_serial[pp](1, 0) * K_serial[pp](2, 2) * N2 * X0 * Circle_params[ii * 6 + 1] - K_serial[pp](1, 2) * K_serial[pp](2, 0) * N2 * X0 * Circle_params[ii * 6 + 1] - K_serial[pp](0, 0) * K_serial[pp](2, 2) * N2 * Y0 * Circle_params[ii * 6 + 1] + K_serial[pp](0, 2) * K_serial[pp](2, 0) * N2 * Y0 * Circle_params[ii * 6 + 1] + K_serial[pp](1, 0) * K_serial[pp](2, 2) * N3 * X0 * Circle_params[ii * 6 + 2] - K_serial[pp](1, 2) * K_serial[pp](2, 0) * N3 * X0 * Circle_params[ii * 6 + 2] - K_serial[pp](0, 0) * K_serial[pp](2, 2) * N3 * Y0 * Circle_params[ii * 6 + 2] + K_serial[pp](0, 2) * K_serial[pp](2, 0) * N3 * Y0 * Circle_params[ii * 6 + 2]) / (K_serial[pp](0, 0) * K_serial[pp](1, 1) * N3 - K_serial[pp](0, 0) * K_serial[pp](1, 2) * N2 - K_serial[pp](0, 1) * K_serial[pp](1, 0) * N3 + K_serial[pp](0, 1) * K_serial[pp](1, 2) * N1 + K_serial[pp](0, 2) * K_serial[pp](1, 0) * N2 - K_serial[pp](0, 2) * K_serial[pp](1, 1) * N1 + K_serial[pp](1, 0) * K_serial[pp](2, 1) * N3 * X0 - K_serial[pp](1, 0) * K_serial[pp](2, 2) * N2 * X0 - K_serial[pp](1, 1) * K_serial[pp](2, 0) * N3 * X0 + K_serial[pp](1, 1) * K_serial[pp](2, 2) * N1 * X0 + K_serial[pp](1, 2) * K_serial[pp](2, 0) * N2 * X0 - K_serial[pp](1, 2) * K_serial[pp](2, 1) * N1 * X0 - K_serial[pp](0, 0) * K_serial[pp](2, 1) * N3 * Y0 + K_serial[pp](0, 0) * K_serial[pp](2, 2) * N2 * Y0 + K_serial[pp](0, 1) * K_serial[pp](2, 0) * N3 * Y0 - K_serial[pp](0, 1) * K_serial[pp](2, 2) * N1 * Y0 - K_serial[pp](0, 2) * K_serial[pp](2, 0) * N2 * Y0 + K_serial[pp](0, 2) * K_serial[pp](2, 1) * N1 * Y0);
						double z_w = (K_serial[pp](0, 0) * K_serial[pp](1, 3) * N2 - K_serial[pp](0, 1) * K_serial[pp](1, 3) * N1 - K_serial[pp](0, 3) * K_serial[pp](1, 0) * N2 + K_serial[pp](0, 3) * K_serial[pp](1, 1) * N1 + K_serial[pp](1, 0) * K_serial[pp](2, 3) * N2 * X0 - K_serial[pp](1, 1) * K_serial[pp](2, 3) * N1 * X0 - K_serial[pp](1, 3) * K_serial[pp](2, 0) * N2 * X0 + K_serial[pp](1, 3) * K_serial[pp](2, 1) * N1 * X0 + K_serial[pp](0, 0) * K_serial[pp](1, 1) * N1 * Circle_params[ii * 6 + 0] - K_serial[pp](0, 1) * K_serial[pp](1, 0) * N1 * Circle_params[ii * 6 + 0] - K_serial[pp](0, 0) * K_serial[pp](2, 3) * N2 * Y0 + K_serial[pp](0, 1) * K_serial[pp](2, 3) * N1 * Y0 + K_serial[pp](0, 3) * K_serial[pp](2, 0) * N2 * Y0 - K_serial[pp](0, 3) * K_serial[pp](2, 1) * N1 * Y0 + K_serial[pp](0, 0) * K_serial[pp](1, 1) * N2 * Circle_params[ii * 6 + 1] - K_serial[pp](0, 1) * K_serial[pp](1, 0) * N2 * Circle_params[ii * 6 + 1] + K_serial[pp](0, 0) * K_serial[pp](1, 1) * N3 * Circle_params[ii * 6 + 2] - K_serial[pp](0, 1) * K_serial[pp](1, 0) * N3 * Circle_params[ii * 6 + 2] + K_serial[pp](1, 0) * K_serial[pp](2, 1) * N1 * X0 * Circle_params[ii * 6 + 0] - K_serial[pp](1, 1) * K_serial[pp](2, 0) * N1 * X0 * Circle_params[ii * 6 + 0] - K_serial[pp](0, 0) * K_serial[pp](2, 1) * N1 * Circle_params[ii * 6 + 0] * Y0 + K_serial[pp](0, 1) * K_serial[pp](2, 0) * N1 * Circle_params[ii * 6 + 0] * Y0 + K_serial[pp](1, 0) * K_serial[pp](2, 1) * N2 * X0 * Circle_params[ii * 6 + 1] - K_serial[pp](1, 1) * K_serial[pp](2, 0) * N2 * X0 * Circle_params[ii * 6 + 1] - K_serial[pp](0, 0) * K_serial[pp](2, 1) * N2 * Y0 * Circle_params[ii * 6 + 1] + K_serial[pp](0, 1) * K_serial[pp](2, 0) * N2 * Y0 * Circle_params[ii * 6 + 1] + K_serial[pp](1, 0) * K_serial[pp](2, 1) * N3 * X0 * Circle_params[ii * 6 + 2] - K_serial[pp](1, 1) * K_serial[pp](2, 0) * N3 * X0 * Circle_params[ii * 6 + 2] - K_serial[pp](0, 0) * K_serial[pp](2, 1) * N3 * Y0 * Circle_params[ii * 6 + 2] + K_serial[pp](0, 1) * K_serial[pp](2, 0) * N3 * Y0 * Circle_params[ii * 6 + 2]) / (K_serial[pp](0, 0) * K_serial[pp](1, 1) * N3 - K_serial[pp](0, 0) * K_serial[pp](1, 2) * N2 - K_serial[pp](0, 1) * K_serial[pp](1, 0) * N3 + K_serial[pp](0, 1) * K_serial[pp](1, 2) * N1 + K_serial[pp](0, 2) * K_serial[pp](1, 0) * N2 - K_serial[pp](0, 2) * K_serial[pp](1, 1) * N1 + K_serial[pp](1, 0) * K_serial[pp](2, 1) * N3 * X0 - K_serial[pp](1, 0) * K_serial[pp](2, 2) * N2 * X0 - K_serial[pp](1, 1) * K_serial[pp](2, 0) * N3 * X0 + K_serial[pp](1, 1) * K_serial[pp](2, 2) * N1 * X0 + K_serial[pp](1, 2) * K_serial[pp](2, 0) * N2 * X0 - K_serial[pp](1, 2) * K_serial[pp](2, 1) * N1 * X0 - K_serial[pp](0, 0) * K_serial[pp](2, 1) * N3 * Y0 + K_serial[pp](0, 0) * K_serial[pp](2, 2) * N2 * Y0 + K_serial[pp](0, 1) * K_serial[pp](2, 0) * N3 * Y0 - K_serial[pp](0, 1) * K_serial[pp](2, 2) * N1 * Y0 - K_serial[pp](0, 2) * K_serial[pp](2, 0) * N2 * Y0 + K_serial[pp](0, 2) * K_serial[pp](2, 1) * N1 * Y0);
						circle_R[ii] += sqrt(pow(x_w - Circle_params[ii * 6 + 0], 2) + pow(y_w - Circle_params[ii * 6 + 1], 2) + pow(z_w - Circle_params[ii * 6 + 2], 2));
						circle_R_Number[ii] += 1.0;
					}
				}
			}
		}
	}
	for (int ii = 0; ii < point_clouds.size(); ii++)
	{
		circle_R[ii] /= circle_R_Number[ii];
		Circle_params[ii * 6 + 3] = circle_R[ii];
	}
	for (int ii = 0; ii < point_clouds.size(); ii++)
	{
		for (int pp = 0; pp < code_circle_serial.size(); pp++)
		{
			for (int qq = 0; qq < code_circle_serial[pp].size(); qq++)
			{
				if (code_circle_serial[pp][qq].code_num == point_clouds[ii].code)
				{
					for (int tt = 0; tt < contour_circle_double_serial[pp][qq].size(); tt++)
					{
						Eigen::Vector2d img(contour_circle_double_serial[pp][qq][tt].x, contour_circle_double_serial[pp][qq][tt].y);
						ceres::CostFunction* cost_function = ProjectErrorCostFunctionCircle_ori::Create(img, K_serial[pp]);
						ceres::LossFunction* loss_function;
						switch (Loss_type)
						{
						case 0:
							loss_function = new ceres::TrivialLoss();
							break;
						case 1:
							loss_function = new ceres::HuberLoss(Loss_value);
							break;
						case 2:
							loss_function = new ceres::SoftLOneLoss(Loss_value);
							break;
						case 3:
							loss_function = new ceres::CauchyLoss(Loss_value);
							break;
						case 4:
							loss_function = new ceres::ArctanLoss(Loss_value);
							break;
						default:
							loss_function = new ceres::CauchyLoss(Loss_value);
							break;
						}
						problem.AddResidualBlock(cost_function, loss_function, &Circle_params[ii * 6]);
					}
				}
			}
		}
	}
	ceres::Solver::Options options;
	options.minimizer_progress_to_stdout = false;
	options.linear_solver_type = ceres::DENSE_NORMAL_CHOLESKY;
	options.max_num_iterations = max_iter_num;
	options.function_tolerance = stop_value;
	options.gradient_tolerance = stop_value;
	options.parameter_tolerance = stop_value;

	ceres::Solve(options, &problem, summary);
	for (int ii = 0; ii < point_clouds.size(); ii++)
	{
		ori_serial[ii].x() = sin(Circle_params[ii * 6 + 4]);
		ori_serial[ii].y() = cos(Circle_params[ii * 6 + 4]) * sin(Circle_params[ii * 6 + 5]);
		ori_serial[ii].z() = cos(Circle_params[ii * 6 + 4]) * cos(Circle_params[ii * 6 + 5]);
		point_clouds[ii].R = Circle_params[ii * 6 + 3];
	}
	delete[] Circle_params;
}
void SFM_optimize::calculate_RT(std::vector<cv::Point3d> point_clouds,
	std::vector<cv::Point2d> img_points,
	Eigen::Quaterniond& camQvec, Eigen::Vector3d& camTvec
	, double* camK, double* camDK, double* camDP, double* camDT
	, unsigned int max_iter_num
	, double stop_value
	, unsigned int num_threads
	, unsigned int timeout
	, unsigned char Loss_type
	, double Loss_value)
{
	Eigen::Matrix<double, 3, 3> I_matrix = Eigen::Matrix<double, 3, 3>::Zero();
	I_matrix(0, 0) = camK[0];
	I_matrix(0, 1) = camK[4];
	I_matrix(0, 2) = camK[2];
	I_matrix(1, 1) = camK[1];
	I_matrix(1, 2) = camK[3];
	I_matrix(2, 2) = 1;

	Eigen::Matrix<double, 1, 6> I_matrix_dis_k = Eigen::Matrix<double, 1, 6>::Zero();
	I_matrix_dis_k(0, 0) = camDK[0];
	I_matrix_dis_k(0, 1) = camDK[1];
	I_matrix_dis_k(0, 2) = camDK[2];
	I_matrix_dis_k(0, 3) = camDK[3];
	I_matrix_dis_k(0, 4) = camDK[4];
	I_matrix_dis_k(0, 5) = camDK[5];


	Eigen::Matrix<double, 1, 2> I_matrix_dis_p = Eigen::Matrix<double, 1, 2>::Zero();
	I_matrix_dis_p(0, 0) = camDP[0];
	I_matrix_dis_p(0, 1) = camDP[1];


	Eigen::Matrix<double, 1, 4> I_matrix_dis_t = Eigen::Matrix<double, 1, 4>::Zero();
	I_matrix_dis_t(0, 0) = camDT[0];
	I_matrix_dis_t(0, 1) = camDT[1];
	I_matrix_dis_t(0, 2) = camDT[2];
	I_matrix_dis_t(0, 3) = camDT[3];
	ceres::Problem problem;
	ceres::LocalParameterization* qvec_parameterization = new ceres::EigenQuaternionParameterization;
	for (int pp = 0; pp < point_clouds.size(); pp++)
	{

		Eigen::Vector2d img(img_points[pp].x, img_points[pp].y);
		Eigen::Vector3d points(point_clouds[pp].x, point_clouds[pp].y, point_clouds[pp].z);
		ceres::CostFunction* cost_function = ProjectErrorCostFunction_PNP::Create(points, img, I_matrix, I_matrix_dis_k, I_matrix_dis_p, I_matrix_dis_t);
		ceres::LossFunction* loss_function;
		switch (Loss_type)
		{
		case 0:
			loss_function = new ceres::TrivialLoss();
			break;
		case 1:
			loss_function = new ceres::HuberLoss(Loss_value);
			break;
		case 2:
			loss_function = new ceres::SoftLOneLoss(Loss_value);
			break;
		case 3:
			loss_function = new ceres::CauchyLoss(Loss_value);
			break;
		case 4:
			loss_function = new ceres::ArctanLoss(Loss_value);
			break;
		default:
			loss_function = new ceres::CauchyLoss(Loss_value);
			break;
		}
		problem.AddResidualBlock(cost_function, loss_function, camQvec.coeffs().data(), camTvec.data());
	}

	problem.SetParameterization(camQvec.coeffs().data(), qvec_parameterization);
	ceres::Solver::Options options;
	options.minimizer_progress_to_stdout = false;
	options.linear_solver_type = ceres::DENSE_SCHUR;
	options.max_num_iterations = max_iter_num;
	options.function_tolerance = stop_value;
	options.gradient_tolerance = stop_value;
	options.parameter_tolerance = stop_value;
	ceres::Solver::Summary summary;
	ceres::Solve(options, &problem, &summary);
}
void SFM_optimize::calculate_SFM(std::vector<std::vector<Coded_detect_inf>> code_circle_serial,
	std::vector<sfm_3d_group> &point_clouds,
	std::vector<Eigen::Quaterniond>& camQvec, std::vector<Eigen::Vector3d>& camTvec
	, std::vector<std::vector<cv::Point2d>>& Re_project_Map
	, double*& camK, double*& camDK, double*& camDP, double*& camDT
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
	, double* C_range)
{
	double* opt_points = new double[3 * point_clouds.size()];
	double* T_sec = new double[2];
	T_sec[0] = asin(camTvec[1][0]);
	T_sec[1] = atan2(camTvec[1][1], camTvec[1][2]);


	Re_project_Map.clear();
	double al_er = 0;
	for (int pp = 0; pp < code_circle_serial.size(); pp++)
	{
		std::vector<cv::Point2d> err_temp;
		for (int qq = 0; qq < code_circle_serial[pp].size(); qq++)
		{
			bool has_exist = false;
			for (int ii = 0; ii < point_clouds.size(); ii++)
			{
				if (code_circle_serial[pp][qq].code_num == point_clouds[ii].code)
				{
					point_clouds[ii].weight++;
					Eigen::Quaternion<double> qvec(camQvec[pp]);
					Eigen::Matrix<double, 3, 1> tvec;
					tvec << camTvec[pp][0], camTvec[pp][1], camTvec[pp][2];

					Eigen::Vector3d obj(point_clouds.at(ii).x, point_clouds.at(ii).y, point_clouds.at(ii).z);
					Eigen::Vector2d img(code_circle_serial[pp][qq].x, code_circle_serial[pp][qq].y);

					Eigen::Matrix<double, 3, 1> obj_cam_coor = qvec.toRotationMatrix() * obj + tvec;

					double a = obj_cam_coor(0) / obj_cam_coor(2);
					double b = obj_cam_coor(1) / obj_cam_coor(2);
					double r2 = (a * a + b * b);
					double r4 = r2 * r2;
					double r6 = r2 * r4;

					double xd = a * (1.0 + camDK[0] * r2 + camDK[1] * r4 + camDK[2] * r6)
						/ (1.0 + camDK[3] * r2 + camDK[4] * r4 + camDK[5] * r6)
						+ 2.0 * camDP[0] * a * b + camDP[1] * (r2 + 2.0 * a * a)
						+ camDT[0] * r2 + camDT[2] * r4;
					double yd = b * (1.0 + camDK[0] * r2 + camDK[1] * r4 + camDK[2] * r6)
						/ (1.0 + camDK[3] * r2 + camDK[4] * r4 + camDK[5] * r6)
						+ 2.0 * camDP[1] * a * b + camDP[0] * (r2 + 2.0 * b * b)
						+ camDT[1] * r2 + camDT[3] * r4;

					double ud = camK[0] * xd + camK[4] * yd + camK[2];
					double vd = camK[1] * yd + camK[3];
					err_temp.push_back(cv::Point2d(ud - img(0), vd - img(1)));

					has_exist = true;
					break;
				}
			}
			if (!has_exist)
			{
				err_temp.push_back(cv::Point2d(std::numeric_limits<double>::max(), std::numeric_limits<double>::max()));
			}
		}
		Re_project_Map.push_back(err_temp);
	}

	for (int ii = 0; ii < point_clouds.size(); ii++)
	{
		point_clouds[ii].weight = 0;
	}
	if (!Use_same_F)
	{
		for (int ii = 0; ii < point_clouds.size(); ii++)
		{
			opt_points[ii * 3] = point_clouds[ii].x;
			opt_points[ii * 3 + 1] = point_clouds[ii].y;
			opt_points[ii * 3 + 2] = point_clouds[ii].z;
		}

		ceres::Problem problem;
		ceres::LocalParameterization* qvec_parameterization = new ceres::EigenQuaternionParameterization;
		for (int ii = 0; ii < point_clouds.size(); ii++)
		{
			for (int pp = 0; pp < code_circle_serial.size(); pp++)
			{
				for (int qq = 0; qq < code_circle_serial[pp].size(); qq++)
				{
					if (code_circle_serial[pp][qq].code_num == point_clouds[ii].code)
					{
						if (pp == 0)
						{
							Eigen::Vector2d img(code_circle_serial[pp][qq].x, code_circle_serial[pp][qq].y);
							ceres::CostFunction* cost_function = ProjectErrorCostFunctionPinehole_SFM_first::Create(img);
							ceres::LossFunction* loss_function;
							switch (Loss_type)
							{
							case 0:
								loss_function = new ceres::TrivialLoss();
								break;
							case 1:
								loss_function = new ceres::HuberLoss(Loss_value);
								break;
							case 2:
								loss_function = new ceres::SoftLOneLoss(Loss_value);
								break;
							case 3:
								loss_function = new ceres::CauchyLoss(Loss_value);
								break;
							case 4:
								loss_function = new ceres::ArctanLoss(Loss_value);
								break;
							default:
								loss_function = new ceres::CauchyLoss(Loss_value);
								break;
							}
							problem.AddResidualBlock(cost_function, loss_function, opt_points + ii * 3, camK, camK + 2, camK + 4,
								camDK, camDK + 1, camDK + 2, camDK + 3, camDK + 4, camDK + 5, camDP, camDT, camDT + 2);
						}
						else if (pp == 1)
						{
							Eigen::Vector2d img(code_circle_serial[pp][qq].x, code_circle_serial[pp][qq].y);
							ceres::CostFunction* cost_function = ProjectErrorCostFunctionPinehole_SFM_Rtheta::Create(img);
							ceres::LossFunction* loss_function;
							switch (Loss_type)
							{
							case 0:
								loss_function = new ceres::TrivialLoss();
								break;
							case 1:
								loss_function = new ceres::HuberLoss(Loss_value);
								break;
							case 2:
								loss_function = new ceres::SoftLOneLoss(Loss_value);
								break;
							case 3:
								loss_function = new ceres::CauchyLoss(Loss_value);
								break;
							case 4:
								loss_function = new ceres::ArctanLoss(Loss_value);
								break;
							default:
								loss_function = new ceres::CauchyLoss(Loss_value);
								break;
							}
							problem.AddResidualBlock(cost_function, loss_function, opt_points + ii * 3, camK, camK + 2, camK + 4,
								camDK, camDK + 1, camDK + 2, camDK + 3, camDK + 4, camDK + 5, camDP, camDT, camDT + 2, camQvec[pp].coeffs().data(), T_sec);
						}
						else
						{
							Eigen::Vector2d img(code_circle_serial[pp][qq].x, code_circle_serial[pp][qq].y);
							ceres::CostFunction* cost_function = ProjectErrorCostFunctionPinehole_SFM::Create(img);
							ceres::LossFunction* loss_function;
							switch (Loss_type)
							{
							case 0:
								loss_function = new ceres::TrivialLoss();
								break;
							case 1:
								loss_function = new ceres::HuberLoss(Loss_value);
								break;
							case 2:
								loss_function = new ceres::SoftLOneLoss(Loss_value);
								break;
							case 3:
								loss_function = new ceres::CauchyLoss(Loss_value);
								break;
							case 4:
								loss_function = new ceres::ArctanLoss(Loss_value);
								break;
							default:
								loss_function = new ceres::CauchyLoss(Loss_value);
								break;
							}
							problem.AddResidualBlock(cost_function, loss_function, opt_points + ii * 3, camK, camK + 2, camK + 4,
								camDK, camDK + 1, camDK + 2, camDK + 3, camDK + 4, camDK + 5, camDP, camDT, camDT + 2, camQvec[pp].coeffs().data(), camTvec[pp].data());
						}
						if (pp != 0)
						{
							problem.SetParameterization(camQvec[pp].coeffs().data(), qvec_parameterization);
						}
					}
				}
			}
		}
		if (Fixed_all)
		{
			problem.SetParameterBlockConstant(camK);
			problem.SetParameterBlockConstant(camK + 2);
			problem.SetParameterBlockConstant(camK + 4);
			for (int tt = 0; tt < 6; tt++)
			{
				problem.SetParameterBlockConstant(camDK + tt);
			}
			problem.SetParameterBlockConstant(camDT);
			problem.SetParameterBlockConstant(camDT + 2);
			problem.SetParameterBlockConstant(camDP);
		}
		else
		{
			if (F_range != nullptr && Use_F)
			{
				if (F_range[2] == 0)
				{
					problem.SetParameterBlockConstant(camK);
				}
				else
				{
					problem.SetParameterLowerBound(camK, 0, F_range[0] - F_range[2]);
					problem.SetParameterLowerBound(camK, 1, F_range[1] - F_range[2]);
					problem.SetParameterUpperBound(camK, 0, F_range[0] + F_range[2]);
					problem.SetParameterUpperBound(camK, 1, F_range[1] + F_range[2]);
				}
			}
			if (C_range != nullptr && Use_Cx_Cy)
			{
				if (C_range[2] == 0)
				{
					problem.SetParameterBlockConstant(camK + 2);
				}
				else
				{
					problem.SetParameterLowerBound(camK + 2, 0, C_range[0] - C_range[2]);
					problem.SetParameterLowerBound(camK + 2, 1, C_range[1] - C_range[2]);
					problem.SetParameterUpperBound(camK + 2, 0, C_range[0] + C_range[2]);
					problem.SetParameterUpperBound(camK + 2, 1, C_range[1] + C_range[2]);
				}
			}
			switch (Dis_P_num)
			{
			case 0:
				camDP[0] = 0;
				camDP[1] = 0;
				problem.SetParameterBlockConstant(camDP);
				break;
			default:
				problem.SetParameterBlockVariable(camDP);
				break;
			}

			switch (Dis_T_num)
			{
			case 0:
				camDT[0] = 0;
				camDT[1] = 0;
				camDT[2] = 0;
				camDT[3] = 0;
				problem.SetParameterBlockConstant(camDT);
				problem.SetParameterBlockConstant(camDT + 2);
				break;
			case 2:
				camDT[2] = 0;
				camDT[3] = 0;
				problem.SetParameterBlockVariable(camDT);
				problem.SetParameterBlockConstant(camDT + 2);
				break;
			default:
				problem.SetParameterBlockVariable(camDT);
				problem.SetParameterBlockVariable(camDT + 2);
				break;
			}

			if (!Use_F)
			{
				problem.SetParameterBlockConstant(camK);
			}
			else
			{
				problem.SetParameterBlockVariable(camK);
			}
			if (!Use_Cx_Cy)
			{
				problem.SetParameterBlockConstant(camK + 2);
			}
			else
			{
				problem.SetParameterBlockVariable(camK + 2);
			}
			if (!Use_shear)
			{
				camK[4] = 0;
				problem.SetParameterBlockConstant(camK + 4);
			}
			else
			{
				problem.SetParameterBlockVariable(camK + 4);
			}
			for (int tt = 0; tt < Dis_K_num; tt++)
			{
				problem.SetParameterBlockVariable(camDK + tt);
			}
			for (int tt = Dis_K_num; tt < 6; tt++)
			{
				camDK[tt] = 0;
				problem.SetParameterBlockConstant(camDK + tt);
			}
		}

		
		ceres::Solver::Options options;
		options.minimizer_progress_to_stdout = false;
		options.linear_solver_type = ceres::DENSE_SCHUR;
		options.max_num_iterations = max_iter_num;
		options.function_tolerance = stop_value;
		options.gradient_tolerance = stop_value;
		options.parameter_tolerance = stop_value;

		ceres::Solve(options, &problem, summary);

		camTvec[1][0] = sin(T_sec[0]);
		camTvec[1][1] = cos(T_sec[0]) * sin(T_sec[1]);
		camTvec[1][2] = cos(T_sec[0]) * cos(T_sec[1]);
		for (int ii = 0; ii < point_clouds.size(); ii++)
		{
			point_clouds[ii].x = opt_points[ii * 3];
			point_clouds[ii].y = opt_points[ii * 3 + 1];
			point_clouds[ii].z = opt_points[ii * 3 + 2];
		}
		Re_project_Map.clear();
		double al_er = 0;
		for (int pp = 0; pp < code_circle_serial.size(); pp++)
		{
			std::vector<cv::Point2d> err_temp;
			for (int qq = 0; qq < code_circle_serial[pp].size(); qq++)
			{
				bool has_exist = false;
				for (int ii = 0; ii < point_clouds.size(); ii++)
				{
					if (code_circle_serial[pp][qq].code_num == point_clouds[ii].code)
					{
						point_clouds[ii].weight++;
						Eigen::Quaternion<double> qvec(camQvec[pp]);
						Eigen::Matrix<double, 3, 1> tvec;
						tvec << camTvec[pp][0], camTvec[pp][1], camTvec[pp][2];

						Eigen::Vector3d obj(point_clouds.at(ii).x, point_clouds.at(ii).y, point_clouds.at(ii).z);
						Eigen::Vector2d img(code_circle_serial[pp][qq].x, code_circle_serial[pp][qq].y);

						Eigen::Matrix<double, 3, 1> obj_cam_coor = qvec.toRotationMatrix() * obj + tvec;

						double a = obj_cam_coor(0) / obj_cam_coor(2);
						double b = obj_cam_coor(1) / obj_cam_coor(2);
						double r2 = (a * a + b * b);
						double r4 = r2 * r2;
						double r6 = r2 * r4;

						double xd = a * (1.0 + camDK[0] * r2 + camDK[1] * r4 + camDK[2] * r6)
							/ (1.0 + camDK[3] * r2 + camDK[4] * r4 + camDK[5] * r6)
							+ 2.0 * camDP[0] * a * b + camDP[1] * (r2 + 2.0 * a * a)
							+ camDT[0] * r2 + camDT[2] * r4;
						double yd = b * (1.0 + camDK[0] * r2 + camDK[1] * r4 + camDK[2] * r6)
							/ (1.0 + camDK[3] * r2 + camDK[4] * r4 + camDK[5] * r6)
							+ 2.0 * camDP[1] * a * b + camDP[0] * (r2 + 2.0 * b * b)
							+ camDT[1] * r2 + camDT[3] * r4;

						double ud = camK[0] * xd + camK[4] * yd + camK[2];
						double vd = camK[1] * yd + camK[3];
						err_temp.push_back(cv::Point2d(ud - img(0), vd - img(1)));
						has_exist = true;
						break;
					}
				}
				if (!has_exist)
				{
					err_temp.push_back(cv::Point2d(std::numeric_limits<double>::max(), std::numeric_limits<double>::max()));
				}
			}
			Re_project_Map.push_back(err_temp);
		}
	}
	else
	{
		for (int ii = 0; ii < point_clouds.size(); ii++)
		{
			opt_points[ii * 3] = point_clouds[ii].x;
			opt_points[ii * 3 + 1] = point_clouds[ii].y;
			opt_points[ii * 3 + 2] = point_clouds[ii].z;
		}

		ceres::Problem problem;
		ceres::LocalParameterization* qvec_parameterization = new ceres::EigenQuaternionParameterization;
		for (int ii = 0; ii < point_clouds.size(); ii++)
		{
			for (int pp = 0; pp < code_circle_serial.size(); pp++)
			{
				for (int qq = 0; qq < code_circle_serial[pp].size(); qq++)
				{
					if (code_circle_serial[pp][qq].code_num == point_clouds[ii].code)
					{
						if (pp == 0)
						{
							Eigen::Vector2d img(code_circle_serial[pp][qq].x, code_circle_serial[pp][qq].y);
							ceres::CostFunction* cost_function = ProjectErrorCostFunctionPinehole_SFM_first_sameF::Create(img);
							ceres::LossFunction* loss_function;
							switch (Loss_type)
							{
							case 0:
								loss_function = new ceres::TrivialLoss();
								break;
							case 1:
								loss_function = new ceres::HuberLoss(Loss_value);
								break;
							case 2:
								loss_function = new ceres::SoftLOneLoss(Loss_value);
								break;
							case 3:
								loss_function = new ceres::CauchyLoss(Loss_value);
								break;
							case 4:
								loss_function = new ceres::ArctanLoss(Loss_value);
								break;
							default:
								loss_function = new ceres::CauchyLoss(Loss_value);
								break;
							}
							problem.AddResidualBlock(cost_function, loss_function, opt_points + ii * 3, camK, camK + 2, camK + 4,
								camDK, camDK + 1, camDK + 2, camDK + 3, camDK + 4, camDK + 5, camDP, camDT, camDT + 2);
						}
						else if (pp == 1)
						{
							Eigen::Vector2d img(code_circle_serial[pp][qq].x, code_circle_serial[pp][qq].y);
							ceres::CostFunction* cost_function = ProjectErrorCostFunctionPinehole_SFM_Rtheta_sameF::Create(img);
							ceres::LossFunction* loss_function;
							switch (Loss_type)
							{
							case 0:
								loss_function = new ceres::TrivialLoss();
								break;
							case 1:
								loss_function = new ceres::HuberLoss(Loss_value);
								break;
							case 2:
								loss_function = new ceres::SoftLOneLoss(Loss_value);
								break;
							case 3:
								loss_function = new ceres::CauchyLoss(Loss_value);
								break;
							case 4:
								loss_function = new ceres::ArctanLoss(Loss_value);
								break;
							default:
								loss_function = new ceres::CauchyLoss(Loss_value);
								break;
							}
							problem.AddResidualBlock(cost_function, loss_function, opt_points + ii * 3, camK, camK + 2, camK + 4,
								camDK, camDK + 1, camDK + 2, camDK + 3, camDK + 4, camDK + 5, camDP, camDT, camDT + 2, camQvec[pp].coeffs().data(), T_sec);
						}
						else
						{
							Eigen::Vector2d img(code_circle_serial[pp][qq].x, code_circle_serial[pp][qq].y);
							ceres::CostFunction* cost_function = ProjectErrorCostFunctionPinehole_SFM_sameF::Create(img);
							ceres::LossFunction* loss_function;
							switch (Loss_type)
							{
							case 0:
								loss_function = new ceres::TrivialLoss();
								break;
							case 1:
								loss_function = new ceres::HuberLoss(Loss_value);
								break;
							case 2:
								loss_function = new ceres::SoftLOneLoss(Loss_value);
								break;
							case 3:
								loss_function = new ceres::CauchyLoss(Loss_value);
								break;
							case 4:
								loss_function = new ceres::ArctanLoss(Loss_value);
								break;
							default:
								loss_function = new ceres::CauchyLoss(Loss_value);
								break;
							}
							problem.AddResidualBlock(cost_function, loss_function, opt_points + ii * 3, camK, camK + 2, camK + 4,
								camDK, camDK + 1, camDK + 2, camDK + 3, camDK + 4, camDK + 5, camDP, camDT, camDT + 2, camQvec[pp].coeffs().data(), camTvec[pp].data());
						}
						if (pp != 0)
						{
							problem.SetParameterization(camQvec[pp].coeffs().data(), qvec_parameterization);
						}
					}
				}
			}
		}


		if (Fixed_all)
		{
			problem.SetParameterBlockConstant(camK);
			problem.SetParameterBlockConstant(camK + 2);
			problem.SetParameterBlockConstant(camK + 4);
			problem.SetParameterBlockConstant(camK);
			problem.SetParameterBlockConstant(camK + 2);
			problem.SetParameterBlockConstant(camK + 4);
			for (int tt = 0; tt < 6; tt++)
			{
				problem.SetParameterBlockConstant(camDK + tt);
			}
			problem.SetParameterBlockConstant(camDT);
			problem.SetParameterBlockConstant(camDT + 2);
			problem.SetParameterBlockConstant(camDP);
		}
		else
		{
			if (F_range != nullptr && Use_F)
			{
				if (F_range[2] == 0)
				{
					problem.SetParameterBlockConstant(camK);
				}
				else
				{
					problem.SetParameterLowerBound(camK, 0, F_range[0] - F_range[2]);
					problem.SetParameterUpperBound(camK, 0, F_range[0] + F_range[2]);
				}
			}
			if (C_range != nullptr && Use_Cx_Cy)
			{
				if (C_range[2] == 0)
				{
					problem.SetParameterBlockConstant(camK + 2);
				}
				else
				{
					problem.SetParameterLowerBound(camK + 2, 0, C_range[0] - C_range[2]);
					problem.SetParameterLowerBound(camK + 2, 1, C_range[1] - C_range[2]);
					problem.SetParameterUpperBound(camK + 2, 0, C_range[0] + C_range[2]);
					problem.SetParameterUpperBound(camK + 2, 1, C_range[1] + C_range[2]);
				}
			}

			switch (Dis_P_num)
			{
			case 0:
				camDP[0] = 0;
				camDP[1] = 0;
				problem.SetParameterBlockConstant(camDP);
				break;
			default:
				problem.SetParameterBlockVariable(camDP);
				break;
			}

			switch (Dis_T_num)
			{
			case 0:
				camDT[0] = 0;
				camDT[1] = 0;
				camDT[2] = 0;
				camDT[3] = 0;
				problem.SetParameterBlockConstant(camDT);
				problem.SetParameterBlockConstant(camDT + 2);
				break;
			case 2:
				camDT[2] = 0;
				camDT[3] = 0;
				problem.SetParameterBlockVariable(camDT);
				problem.SetParameterBlockConstant(camDT + 2);
				break;
			default:
				problem.SetParameterBlockVariable(camDT);
				problem.SetParameterBlockVariable(camDT + 2);
				break;
			}

			if (!Use_F)
			{
				problem.SetParameterBlockConstant(camK);
			}
			else
			{
				problem.SetParameterBlockVariable(camK);
			}
			if (!Use_Cx_Cy)
			{
				problem.SetParameterBlockConstant(camK + 2);
			}
			else
			{
				problem.SetParameterBlockVariable(camK + 2);
			}
			if (!Use_shear)
			{
				camK[4] = 0;
				problem.SetParameterBlockConstant(camK + 4);
			}
			else
			{
				problem.SetParameterBlockVariable(camK + 4);
			}
			for (int tt = 0; tt < Dis_K_num; tt++)
			{
				problem.SetParameterBlockVariable(camDK + tt);
			}
			for (int tt = Dis_K_num; tt < 6; tt++)
			{
				camDK[tt] = 0;
				problem.SetParameterBlockConstant(camDK + tt);
			}
		}
		ceres::Solver::Options options;
		options.minimizer_progress_to_stdout = false;
		options.linear_solver_type = ceres::DENSE_SCHUR;
		options.max_num_iterations = max_iter_num;
		options.function_tolerance = stop_value;
		options.gradient_tolerance = stop_value;
		options.parameter_tolerance = stop_value;
		ceres::Solve(options, &problem, summary);

		camK[1] = camK[0];
		camTvec[1][0] = sin(T_sec[0]);
		camTvec[1][1] = cos(T_sec[0]) * sin(T_sec[1]);
		camTvec[1][2] = cos(T_sec[0]) * cos(T_sec[1]);
		for (int ii = 0; ii < point_clouds.size(); ii++)
		{
			point_clouds[ii].x = opt_points[ii * 3];
			point_clouds[ii].y = opt_points[ii * 3 + 1];
			point_clouds[ii].z = opt_points[ii * 3 + 2];
		}

		Re_project_Map.clear();
		double al_er = 0;
		for (int pp = 0; pp < code_circle_serial.size(); pp++)
		{
			std::vector<cv::Point2d> err_temp;
			for (int qq = 0; qq < code_circle_serial[pp].size(); qq++)
			{
				bool has_exist = false;
				for (int ii = 0; ii < point_clouds.size(); ii++)
				{
					if (code_circle_serial[pp][qq].code_num == point_clouds[ii].code)
					{
						point_clouds[ii].weight++;
						Eigen::Quaternion<double> qvec(camQvec[pp]);
						Eigen::Matrix<double, 3, 1> tvec;
						tvec << camTvec[pp][0], camTvec[pp][1], camTvec[pp][2];

						Eigen::Vector3d obj(point_clouds.at(ii).x, point_clouds.at(ii).y, point_clouds.at(ii).z);
						Eigen::Vector2d img(code_circle_serial[pp][qq].x, code_circle_serial[pp][qq].y);

						Eigen::Matrix<double, 3, 1> obj_cam_coor = qvec.toRotationMatrix() * obj + tvec;

						double a = obj_cam_coor(0) / obj_cam_coor(2);
						double b = obj_cam_coor(1) / obj_cam_coor(2);
						double r2 = (a * a + b * b);
						double r4 = r2 * r2;
						double r6 = r2 * r4;

						double xd = a * (1.0 + camDK[0] * r2 + camDK[1] * r4 + camDK[2] * r6)
							/ (1.0 + camDK[3] * r2 + camDK[4] * r4 + camDK[5] * r6)
							+ 2.0 * camDP[0] * a * b + camDP[1] * (r2 + 2.0 * a * a)
							+ camDT[0] * r2 + camDT[2] * r4;
						double yd = b * (1.0 + camDK[0] * r2 + camDK[1] * r4 + camDK[2] * r6)
							/ (1.0 + camDK[3] * r2 + camDK[4] * r4 + camDK[5] * r6)
							+ 2.0 * camDP[1] * a * b + camDP[0] * (r2 + 2.0 * b * b)
							+ camDT[1] * r2 + camDT[3] * r4;

						double ud = camK[0] * xd + camK[4] * yd + camK[2];
						double vd = camK[1] * yd + camK[3];
						err_temp.push_back(cv::Point2d(ud - img(0), vd - img(1)));

						has_exist = true;
						break;
					}
				}
				if (!has_exist)
				{
					err_temp.push_back(cv::Point2d(std::numeric_limits<double>::max(), std::numeric_limits<double>::max()));
				}
			}
			Re_project_Map.push_back(err_temp);
		}
	}
}


void SFM_optimize::calculate_SFM_part(std::vector<std::vector<Coded_detect_inf>> code_circle_serial,
	std::vector<sfm_3d_group>& point_clouds,
	std::vector<Eigen::Quaterniond>& camQvec, std::vector<Eigen::Vector3d>& camTvec
	, std::vector<std::vector<cv::Point2d>>& Re_project_Map
	, double*& camK, double*& camDK, double*& camDP, double*& camDT
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
	, double* C_range)
{
	double* opt_points = new double[3 * point_clouds.size()];
	double* T_sec = new double[2];
	T_sec[0] = asin(camTvec[1][0]);
	T_sec[1] = atan2(camTvec[1][1], camTvec[1][2]);

	Re_project_Map.clear();
	double al_er = 0;
	for (int pp = 0; pp < code_circle_serial.size(); pp++)
	{
		std::vector<cv::Point2d> err_temp;
		for (int qq = 0; qq < code_circle_serial[pp].size(); qq++)
		{
			bool has_exist = false;
			for (int ii = 0; ii < point_clouds.size(); ii++)
			{
				if (code_circle_serial[pp][qq].code_num == point_clouds[ii].code)
				{
					point_clouds[ii].weight++;
					Eigen::Quaternion<double> qvec(camQvec[pp]);
					Eigen::Matrix<double, 3, 1> tvec;
					tvec << camTvec[pp][0], camTvec[pp][1], camTvec[pp][2];

					Eigen::Vector3d obj(point_clouds.at(ii).x, point_clouds.at(ii).y, point_clouds.at(ii).z);
					Eigen::Vector2d img(code_circle_serial[pp][qq].x, code_circle_serial[pp][qq].y);

					Eigen::Matrix<double, 3, 1> obj_cam_coor = qvec.toRotationMatrix() * obj + tvec;

					double a = obj_cam_coor(0) / obj_cam_coor(2);
					double b = obj_cam_coor(1) / obj_cam_coor(2);
					double r2 = (a * a + b * b);
					double r4 = r2 * r2;
					double r6 = r2 * r4;

					double xd = a * (1.0 + camDK[0] * r2 + camDK[1] * r4 + camDK[2] * r6)
						/ (1.0 + camDK[3] * r2 + camDK[4] * r4 + camDK[5] * r6)
						+ 2.0 * camDP[0] * a * b + camDP[1] * (r2 + 2.0 * a * a)
						+ camDT[0] * r2 + camDT[2] * r4;
					double yd = b * (1.0 + camDK[0] * r2 + camDK[1] * r4 + camDK[2] * r6)
						/ (1.0 + camDK[3] * r2 + camDK[4] * r4 + camDK[5] * r6)
						+ 2.0 * camDP[1] * a * b + camDP[0] * (r2 + 2.0 * b * b)
						+ camDT[1] * r2 + camDT[3] * r4;

					double ud = camK[0] * xd + camK[4] * yd + camK[2];
					double vd = camK[1] * yd + camK[3];
					err_temp.push_back(cv::Point2d(ud - img(0), vd - img(1)));
					has_exist = true;
					break;
				}
			}
			if (!has_exist)
			{
				err_temp.push_back(cv::Point2d(std::numeric_limits<double>::max(), std::numeric_limits<double>::max()));
			}
		}
		Re_project_Map.push_back(err_temp);
	}


	std::vector<bool> poi_enable;
	std::vector<bool> view_enable;
	for (int ii = 0; ii < point_clouds.size(); ii++)
	{
		bool is_ena = false;
		for (int jj = 0; jj < code_circle_serial[code_circle_serial.size() - 1].size(); jj++)
		{
			if (point_clouds[ii].code == code_circle_serial[code_circle_serial.size() - 1][jj].code_num)
			{
				is_ena = true;
				break;
			}
		}
		poi_enable.push_back(is_ena);
	}

	for (int ii = 0; ii < code_circle_serial.size() - 1; ii++)
	{
		view_enable.push_back(false);
	}
	view_enable.push_back(true);
	for (int ii = 0; ii < point_clouds.size(); ii++)
	{
		point_clouds[ii].weight = 0;
	}
	if (!Use_same_F)
	{
		for (int ii = 0; ii < point_clouds.size(); ii++)
		{
			opt_points[ii * 3] = point_clouds[ii].x;
			opt_points[ii * 3 + 1] = point_clouds[ii].y;
			opt_points[ii * 3 + 2] = point_clouds[ii].z;
		}

		ceres::Problem problem;
		ceres::LocalParameterization* qvec_parameterization = new ceres::EigenQuaternionParameterization;
		for (int ii = 0; ii < point_clouds.size(); ii++)
		{
			for (int pp = 0; pp < code_circle_serial.size(); pp++)
			{
				for (int qq = 0; qq < code_circle_serial[pp].size(); qq++)
				{
					if (code_circle_serial[pp][qq].code_num == point_clouds[ii].code)
					{
						if (pp == 0)
						{
							Eigen::Vector2d img(code_circle_serial[pp][qq].x, code_circle_serial[pp][qq].y);
							ceres::CostFunction* cost_function = ProjectErrorCostFunctionPinehole_SFM_first::Create(img);
							ceres::LossFunction* loss_function;
							switch (Loss_type)
							{
							case 0:
								loss_function = new ceres::TrivialLoss();
								break;
							case 1:
								loss_function = new ceres::HuberLoss(Loss_value);
								break;
							case 2:
								loss_function = new ceres::SoftLOneLoss(Loss_value);
								break;
							case 3:
								loss_function = new ceres::CauchyLoss(Loss_value);
								break;
							case 4:
								loss_function = new ceres::ArctanLoss(Loss_value);
								break;
							default:
								loss_function = new ceres::CauchyLoss(Loss_value);
								break;
							}
							problem.AddResidualBlock(cost_function, loss_function, opt_points + ii * 3, camK, camK + 2, camK + 4,
								camDK, camDK + 1, camDK + 2, camDK + 3, camDK + 4, camDK + 5, camDP, camDT, camDT + 2);
						}
						else if (pp == 1)
						{
							Eigen::Vector2d img(code_circle_serial[pp][qq].x, code_circle_serial[pp][qq].y);
							ceres::CostFunction* cost_function = ProjectErrorCostFunctionPinehole_SFM_Rtheta::Create(img);
							ceres::LossFunction* loss_function;
							switch (Loss_type)
							{
							case 0:
								loss_function = new ceres::TrivialLoss();
								break;
							case 1:
								loss_function = new ceres::HuberLoss(Loss_value);
								break;
							case 2:
								loss_function = new ceres::SoftLOneLoss(Loss_value);
								break;
							case 3:
								loss_function = new ceres::CauchyLoss(Loss_value);
								break;
							case 4:
								loss_function = new ceres::ArctanLoss(Loss_value);
								break;
							default:
								loss_function = new ceres::CauchyLoss(Loss_value);
								break;
							}
							problem.AddResidualBlock(cost_function, loss_function, opt_points + ii * 3, camK, camK + 2, camK + 4,
								camDK, camDK + 1, camDK + 2, camDK + 3, camDK + 4, camDK + 5, camDP, camDT, camDT + 2, camQvec[pp].coeffs().data(), T_sec);
						}
						else
						{
							Eigen::Vector2d img(code_circle_serial[pp][qq].x, code_circle_serial[pp][qq].y);
							ceres::CostFunction* cost_function = ProjectErrorCostFunctionPinehole_SFM::Create(img);
							ceres::LossFunction* loss_function;
							switch (Loss_type)
							{
							case 0:
								loss_function = new ceres::TrivialLoss();
								break;
							case 1:
								loss_function = new ceres::HuberLoss(Loss_value);
								break;
							case 2:
								loss_function = new ceres::SoftLOneLoss(Loss_value);
								break;
							case 3:
								loss_function = new ceres::CauchyLoss(Loss_value);
								break;
							case 4:
								loss_function = new ceres::ArctanLoss(Loss_value);
								break;
							default:
								loss_function = new ceres::CauchyLoss(Loss_value);
								break;
							}
							problem.AddResidualBlock(cost_function, loss_function, opt_points + ii * 3, camK, camK + 2, camK + 4,
								camDK, camDK + 1, camDK + 2, camDK + 3, camDK + 4, camDK + 5, camDP, camDT, camDT + 2, camQvec[pp].coeffs().data(), camTvec[pp].data());
						}
						if (pp != 0)
						{
							problem.SetParameterization(camQvec[pp].coeffs().data(), qvec_parameterization);
						}
					}
				}
			}
		}
		problem.SetParameterBlockConstant(camK);
		problem.SetParameterBlockConstant(camK + 2);
		problem.SetParameterBlockConstant(camK + 4);

		for (int tt = 0; tt < 6; tt++)
		{
			problem.SetParameterBlockConstant(camDK + tt);
		}
		problem.SetParameterBlockConstant(camDT);
		problem.SetParameterBlockConstant(camDT + 2);
		problem.SetParameterBlockConstant(camDP);

		for (int tt = 0; tt < point_clouds.size(); tt++)
		{
			if (!poi_enable[tt])
			{
				problem.SetParameterBlockConstant(opt_points + tt * 3);
			}
		}

		for (int tt = 0; tt < code_circle_serial.size() - 1; tt++)
		{
			if (!view_enable[tt])
			{
				if (tt == 0)
				{

				}
				else if (tt == 1)
				{
					problem.SetParameterBlockConstant(camQvec[tt].coeffs().data());
					problem.SetParameterBlockConstant(T_sec);
				}
				else
				{
					problem.SetParameterBlockConstant(camQvec[tt].coeffs().data());
					problem.SetParameterBlockConstant(camTvec[tt].data());
				}
			}
		}

		ceres::Solver::Options options;
		options.minimizer_progress_to_stdout = false;
		options.linear_solver_type = ceres::DENSE_SCHUR;
		options.max_num_iterations = max_iter_num;
		options.function_tolerance = stop_value;
		options.gradient_tolerance = stop_value;
		options.parameter_tolerance = stop_value;

		ceres::Solve(options, &problem, summary);

		camTvec[1][0] = sin(T_sec[0]);
		camTvec[1][1] = cos(T_sec[0]) * sin(T_sec[1]);
		camTvec[1][2] = cos(T_sec[0]) * cos(T_sec[1]);
		for (int ii = 0; ii < point_clouds.size(); ii++)
		{
			point_clouds[ii].x = opt_points[ii * 3];
			point_clouds[ii].y = opt_points[ii * 3 + 1];
			point_clouds[ii].z = opt_points[ii * 3 + 2];
		}

		Re_project_Map.clear();
		double al_er = 0;
		for (int pp = 0; pp < code_circle_serial.size(); pp++)
		{
			std::vector<cv::Point2d> err_temp;
			for (int qq = 0; qq < code_circle_serial[pp].size(); qq++)
			{
				bool has_exist = false;
				for (int ii = 0; ii < point_clouds.size(); ii++)
				{
					if (code_circle_serial[pp][qq].code_num == point_clouds[ii].code)
					{
						point_clouds[ii].weight++;
						Eigen::Quaternion<double> qvec(camQvec[pp]);
						Eigen::Matrix<double, 3, 1> tvec;
						tvec << camTvec[pp][0], camTvec[pp][1], camTvec[pp][2];

						Eigen::Vector3d obj(point_clouds.at(ii).x, point_clouds.at(ii).y, point_clouds.at(ii).z);
						Eigen::Vector2d img(code_circle_serial[pp][qq].x, code_circle_serial[pp][qq].y);

						Eigen::Matrix<double, 3, 1> obj_cam_coor = qvec.toRotationMatrix() * obj + tvec;

						double a = obj_cam_coor(0) / obj_cam_coor(2);
						double b = obj_cam_coor(1) / obj_cam_coor(2);
						double r2 = (a * a + b * b);
						double r4 = r2 * r2;
						double r6 = r2 * r4;

						double xd = a * (1.0 + camDK[0] * r2 + camDK[1] * r4 + camDK[2] * r6)
							/ (1.0 + camDK[3] * r2 + camDK[4] * r4 + camDK[5] * r6)
							+ 2.0 * camDP[0] * a * b + camDP[1] * (r2 + 2.0 * a * a)
							+ camDT[0] * r2 + camDT[2] * r4;
						double yd = b * (1.0 + camDK[0] * r2 + camDK[1] * r4 + camDK[2] * r6)
							/ (1.0 + camDK[3] * r2 + camDK[4] * r4 + camDK[5] * r6)
							+ 2.0 * camDP[1] * a * b + camDP[0] * (r2 + 2.0 * b * b)
							+ camDT[1] * r2 + camDT[3] * r4;

						double ud = camK[0] * xd + camK[4] * yd + camK[2];
						double vd = camK[1] * yd + camK[3];
						err_temp.push_back(cv::Point2d(ud - img(0), vd - img(1)));
						has_exist = true;
						break;
					}
				}
				if (!has_exist)
				{
					err_temp.push_back(cv::Point2d(std::numeric_limits<double>::max(), std::numeric_limits<double>::max()));
				}
			}
			Re_project_Map.push_back(err_temp);
		}
	}
	else
	{
		for (int ii = 0; ii < point_clouds.size(); ii++)
		{
			opt_points[ii * 3] = point_clouds[ii].x;
			opt_points[ii * 3 + 1] = point_clouds[ii].y;
			opt_points[ii * 3 + 2] = point_clouds[ii].z;
		}

		ceres::Problem problem;
		ceres::LocalParameterization* qvec_parameterization = new ceres::EigenQuaternionParameterization;
		for (int ii = 0; ii < point_clouds.size(); ii++)
		{
			for (int pp = 0; pp < code_circle_serial.size(); pp++)
			{
				for (int qq = 0; qq < code_circle_serial[pp].size(); qq++)
				{
					if (code_circle_serial[pp][qq].code_num == point_clouds[ii].code)
					{
						if (pp == 0)
						{
							Eigen::Vector2d img(code_circle_serial[pp][qq].x, code_circle_serial[pp][qq].y);
							ceres::CostFunction* cost_function = ProjectErrorCostFunctionPinehole_SFM_first_sameF::Create(img);
							ceres::LossFunction* loss_function;
							switch (Loss_type)
							{
							case 0:
								loss_function = new ceres::TrivialLoss();
								break;
							case 1:
								loss_function = new ceres::HuberLoss(Loss_value);
								break;
							case 2:
								loss_function = new ceres::SoftLOneLoss(Loss_value);
								break;
							case 3:
								loss_function = new ceres::CauchyLoss(Loss_value);
								break;
							case 4:
								loss_function = new ceres::ArctanLoss(Loss_value);
								break;
							default:
								loss_function = new ceres::CauchyLoss(Loss_value);
								break;
							}
							problem.AddResidualBlock(cost_function, loss_function, opt_points + ii * 3, camK, camK + 2, camK + 4,
								camDK, camDK + 1, camDK + 2, camDK + 3, camDK + 4, camDK + 5, camDP, camDT, camDT + 2);
						}
						else if (pp == 1)
						{
							Eigen::Vector2d img(code_circle_serial[pp][qq].x, code_circle_serial[pp][qq].y);
							ceres::CostFunction* cost_function = ProjectErrorCostFunctionPinehole_SFM_Rtheta_sameF::Create(img);
							ceres::LossFunction* loss_function;
							switch (Loss_type)
							{
							case 0:
								loss_function = new ceres::TrivialLoss();
								break;
							case 1:
								loss_function = new ceres::HuberLoss(Loss_value);
								break;
							case 2:
								loss_function = new ceres::SoftLOneLoss(Loss_value);
								break;
							case 3:
								loss_function = new ceres::CauchyLoss(Loss_value);
								break;
							case 4:
								loss_function = new ceres::ArctanLoss(Loss_value);
								break;
							default:
								loss_function = new ceres::CauchyLoss(Loss_value);
								break;
							}
							problem.AddResidualBlock(cost_function, loss_function, opt_points + ii * 3, camK, camK + 2, camK + 4,
								camDK, camDK + 1, camDK + 2, camDK + 3, camDK + 4, camDK + 5, camDP, camDT, camDT + 2, camQvec[pp].coeffs().data(), T_sec);
						}
						else
						{
							Eigen::Vector2d img(code_circle_serial[pp][qq].x, code_circle_serial[pp][qq].y);
							ceres::CostFunction* cost_function = ProjectErrorCostFunctionPinehole_SFM_sameF::Create(img);
							ceres::LossFunction* loss_function;
							switch (Loss_type)
							{
							case 0:
								loss_function = new ceres::TrivialLoss();
								break;
							case 1:
								loss_function = new ceres::HuberLoss(Loss_value);
								break;
							case 2:
								loss_function = new ceres::SoftLOneLoss(Loss_value);
								break;
							case 3:
								loss_function = new ceres::CauchyLoss(Loss_value);
								break;
							case 4:
								loss_function = new ceres::ArctanLoss(Loss_value);
								break;
							default:
								loss_function = new ceres::CauchyLoss(Loss_value);
								break;
							}
							problem.AddResidualBlock(cost_function, loss_function, opt_points + ii * 3, camK, camK + 2, camK + 4,
								camDK, camDK + 1, camDK + 2, camDK + 3, camDK + 4, camDK + 5, camDP, camDT, camDT + 2, camQvec[pp].coeffs().data(), camTvec[pp].data());
						}
						if (pp != 0)
						{
							problem.SetParameterization(camQvec[pp].coeffs().data(), qvec_parameterization);
						}
					}
				}
			}
		}


		problem.SetParameterBlockConstant(camK);
		problem.SetParameterBlockConstant(camK + 2);
		problem.SetParameterBlockConstant(camK + 4);
		problem.SetParameterBlockConstant(camK);
		problem.SetParameterBlockConstant(camK + 2);
		problem.SetParameterBlockConstant(camK + 4);
		for (int tt = 0; tt < 6; tt++)
		{
			problem.SetParameterBlockConstant(camDK + tt);
		}
		problem.SetParameterBlockConstant(camDT);
		problem.SetParameterBlockConstant(camDT + 2);
		problem.SetParameterBlockConstant(camDP);

		for (int tt = 0; tt < point_clouds.size(); tt++)
		{
			if (!poi_enable[tt])
			{
				problem.SetParameterBlockConstant(opt_points + tt * 3);
			}
		}

		for (int tt = 0; tt < code_circle_serial.size(); tt++)
		{
			if (!view_enable[tt])
			{
				if (tt == 0)
				{

				}
				else if (tt == 1)
				{
					problem.SetParameterBlockConstant(camQvec[tt].coeffs().data());
					problem.SetParameterBlockConstant(T_sec);
				}
				else
				{
					problem.SetParameterBlockConstant(camQvec[tt].coeffs().data());
					problem.SetParameterBlockConstant(camTvec[tt].data());
				}
			}
		}

		ceres::Solver::Options options;
		options.minimizer_progress_to_stdout = false;
		options.linear_solver_type = ceres::DENSE_SCHUR;
		options.max_num_iterations = max_iter_num;
		options.function_tolerance = stop_value;
		options.gradient_tolerance = stop_value;
		options.parameter_tolerance = stop_value;
		ceres::Solve(options, &problem, summary);

		camK[1] = camK[0];
		camTvec[1][0] = sin(T_sec[0]);
		camTvec[1][1] = cos(T_sec[0]) * sin(T_sec[1]);
		camTvec[1][2] = cos(T_sec[0]) * cos(T_sec[1]);
		for (int ii = 0; ii < point_clouds.size(); ii++)
		{
			point_clouds[ii].x = opt_points[ii * 3];
			point_clouds[ii].y = opt_points[ii * 3 + 1];
			point_clouds[ii].z = opt_points[ii * 3 + 2];
		}

		Re_project_Map.clear();
		double al_er = 0;
		for (int pp = 0; pp < code_circle_serial.size(); pp++)
		{
			std::vector<cv::Point2d> err_temp;
			for (int qq = 0; qq < code_circle_serial[pp].size(); qq++)
			{
				bool has_exist = false;
				for (int ii = 0; ii < point_clouds.size(); ii++)
				{
					if (code_circle_serial[pp][qq].code_num == point_clouds[ii].code)
					{
						point_clouds[ii].weight++;
						Eigen::Quaternion<double> qvec(camQvec[pp]);
						Eigen::Matrix<double, 3, 1> tvec;
						tvec << camTvec[pp][0], camTvec[pp][1], camTvec[pp][2];

						Eigen::Vector3d obj(point_clouds.at(ii).x, point_clouds.at(ii).y, point_clouds.at(ii).z);
						Eigen::Vector2d img(code_circle_serial[pp][qq].x, code_circle_serial[pp][qq].y);

						Eigen::Matrix<double, 3, 1> obj_cam_coor = qvec.toRotationMatrix() * obj + tvec;

						double a = obj_cam_coor(0) / obj_cam_coor(2);
						double b = obj_cam_coor(1) / obj_cam_coor(2);
						double r2 = (a * a + b * b);
						double r4 = r2 * r2;
						double r6 = r2 * r4;

						double xd = a * (1.0 + camDK[0] * r2 + camDK[1] * r4 + camDK[2] * r6)
							/ (1.0 + camDK[3] * r2 + camDK[4] * r4 + camDK[5] * r6)
							+ 2.0 * camDP[0] * a * b + camDP[1] * (r2 + 2.0 * a * a)
							+ camDT[0] * r2 + camDT[2] * r4;
						double yd = b * (1.0 + camDK[0] * r2 + camDK[1] * r4 + camDK[2] * r6)
							/ (1.0 + camDK[3] * r2 + camDK[4] * r4 + camDK[5] * r6)
							+ 2.0 * camDP[1] * a * b + camDP[0] * (r2 + 2.0 * b * b)
							+ camDT[1] * r2 + camDT[3] * r4;

						double ud = camK[0] * xd + camK[4] * yd + camK[2];
						double vd = camK[1] * yd + camK[3];
						err_temp.push_back(cv::Point2d(ud - img(0), vd - img(1)));
						has_exist = true;
						break;
					}
				}
				if (!has_exist)
				{
					err_temp.push_back(cv::Point2d(std::numeric_limits<double>::max(), std::numeric_limits<double>::max()));
				}
			}
			Re_project_Map.push_back(err_temp);
		}
	}
}

void SFM_optimize::calculate_dual(std::vector<std::vector<Coded_detect_inf>> code_circle_serial,
	std::vector<sfm_3d_group>& point_clouds,
	std::vector<Eigen::Quaterniond>& camQvec, std::vector<Eigen::Vector3d>& camTvec
	, std::vector<std::vector<cv::Point2d>>& Re_project_Map
	, double*& camK, double*& camDK, double*& camDP, double*& camDT
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
	, double* C_range)
{
	double* opt_points = new double[3 * point_clouds.size()];
	double* T_sec = new double[2];
	T_sec[0] = asin(camTvec[1][0]);
	T_sec[1] = atan2(camTvec[1][1], camTvec[1][2]);


	cv::Mat K(cv::Matx33d(
		camK[0], 0, camK[2],
		0, camK[1], camK[3],
		0, 0, 1));
	cv::Mat mat_temp_L;
	cv::eigen2cv(camQvec[0].matrix(), mat_temp_L);
	cv::Mat mat_temp2_L;
	cv::eigen2cv(camTvec[0], mat_temp2_L);
	cv::Mat mat_temp_R;
	cv::eigen2cv(camQvec[1].matrix(), mat_temp_R);
	cv::Mat mat_temp2_R;
	cv::eigen2cv(camTvec[1], mat_temp2_R);

	cv::Mat proj1(3, 4, CV_64FC1);
	cv::Mat proj2(3, 4, CV_64FC1);
	for (int ii = 0; ii < 3; ii++)
	{
		for (int jj = 0; jj < 3; jj++)
		{
			proj1.at<double>(ii, jj) = mat_temp_L.at<double>(ii, jj);
			proj2.at<double>(ii, jj) = mat_temp_R.at<double>(ii, jj);
		}
		proj1.at<double>(ii, 3) = mat_temp2_L.at<double>(ii, 0);
		proj2.at<double>(ii, 3) = mat_temp2_R.at<double>(ii, 0);
	}

	proj1 = K * proj1;
	proj2 = K * proj2;
	std::vector<cv::Point2d> P1;
	std::vector<cv::Point2d> P2;
	std::vector<int> Code;
	for (int ii = 0; ii < code_circle_serial[0].size(); ii++)
	{
		for (int jj = 0; jj < code_circle_serial[1].size(); jj++)
		{
			if (code_circle_serial[1][jj].code_num == code_circle_serial[0][ii].code_num)
			{
				P1.push_back(cv::Point2d(code_circle_serial[0][ii].x, code_circle_serial[0][ii].y));
				P2.push_back(cv::Point2d(code_circle_serial[1][jj].x, code_circle_serial[1][jj].y));
				Code.push_back(code_circle_serial[1][jj].code_num);
			}
		}
	}
	std::vector<sfm_3d_group> d3_group;
	cv::Mat structure;
	triangulatePoints(proj1, proj2, P1, P2, structure);
	point_clouds.clear();
	for (int ii = 0; ii < P1.size(); ii++)
	{
		double x = structure.at<double>(0, ii) / structure.at<double>(3, ii);
		double y = structure.at<double>(1, ii) / structure.at<double>(3, ii);
		double z = structure.at<double>(2, ii) / structure.at<double>(3, ii);
		point_clouds.push_back(sfm_3d_group(Code[ii], 100, x, y, z, 0));
	}

	Re_project_Map.clear();
	double al_er = 0;
	for (int pp = 0; pp < code_circle_serial.size(); pp++)
	{
		std::vector<cv::Point2d> err_temp;
		for (int qq = 0; qq < code_circle_serial[pp].size(); qq++)
		{
			bool has_exist = false;
			for (int ii = 0; ii < point_clouds.size(); ii++)
			{
				if (code_circle_serial[pp][qq].code_num == point_clouds[ii].code)
				{
					point_clouds[ii].weight++;
					Eigen::Quaternion<double> qvec(camQvec[pp]);
					Eigen::Matrix<double, 3, 1> tvec;
					tvec << camTvec[pp][0], camTvec[pp][1], camTvec[pp][2];

					Eigen::Vector3d obj(point_clouds.at(ii).x, point_clouds.at(ii).y, point_clouds.at(ii).z);
					Eigen::Vector2d img(code_circle_serial[pp][qq].x, code_circle_serial[pp][qq].y);

					Eigen::Matrix<double, 3, 1> obj_cam_coor = qvec.toRotationMatrix() * obj + tvec;

					double a = obj_cam_coor(0) / obj_cam_coor(2);
					double b = obj_cam_coor(1) / obj_cam_coor(2);
					double r2 = (a * a + b * b);
					double r4 = r2 * r2;
					double r6 = r2 * r4;

					double xd = a * (1.0 + camDK[0] * r2 + camDK[1] * r4 + camDK[2] * r6)
						/ (1.0 + camDK[3] * r2 + camDK[4] * r4 + camDK[5] * r6)
						+ 2.0 * camDP[0] * a * b + camDP[1] * (r2 + 2.0 * a * a)
						+ camDT[0] * r2 + camDT[2] * r4;
					double yd = b * (1.0 + camDK[0] * r2 + camDK[1] * r4 + camDK[2] * r6)
						/ (1.0 + camDK[3] * r2 + camDK[4] * r4 + camDK[5] * r6)
						+ 2.0 * camDP[1] * a * b + camDP[0] * (r2 + 2.0 * b * b)
						+ camDT[1] * r2 + camDT[3] * r4;

					double ud = camK[0] * xd + camK[4] * yd + camK[2];
					double vd = camK[1] * yd + camK[3];
					err_temp.push_back(cv::Point2d(ud - img(0), vd - img(1)));
					has_exist = true;
					break;
				}
			}
			if (!has_exist)
			{
				err_temp.push_back(cv::Point2d(std::numeric_limits<double>::max(), std::numeric_limits<double>::max()));
			}
		}
		Re_project_Map.push_back(err_temp);
	}
}