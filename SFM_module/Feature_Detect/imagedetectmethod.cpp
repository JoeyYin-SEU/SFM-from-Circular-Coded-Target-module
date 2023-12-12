#include "imagedetectmethod.h"
#include <imgproc/types_c.h>

ImageDetectMethod::ImageDetectMethod(QObject* parent)
	: QObject(parent)
{
}

ImageDetectMethod::~ImageDetectMethod()
{
}


bool ImageDetectMethod::ImagePreprocess(const cv::Mat& ori_image_mat, cv::Mat& processed_image_mat)
{
	/////////高斯滤波
	//cv::bilateralFilter(ori_image_mat, processed_image_mat, 9, 50, 25 / 2);
	GaussianBlur(ori_image_mat, processed_image_mat, cv::Size(5, 5), 1, 1);
	return true;
}
bool ImageDetectMethod::DetectClosedContours(const cv::Mat& ori_image_mat, std::vector<std::vector<cv::Point>>& contours, DetectContoursMethod image_process_method)
{
	//vector<vector<cv::Point>> contours;
	std::vector<cv::Vec4i> hierarchy;

	cv::Mat threshold_image;
	if (image_process_method == Sobel_Method)
	{
		threshold_image = ori_image_mat.clone();
		if (threshold_image.channels() == 3)
		{
			cv::cvtColor(threshold_image, threshold_image, CV_BGR2GRAY);
		}
		coder::array<unsigned char, 2U> binary_I;
		binary_I.set_size(threshold_image.rows, threshold_image.cols);
		coder::array<boolean_T, 2U> b_I;
		b_I.set_size(threshold_image.rows, threshold_image.cols);
		for (int idx0{ 0 }; idx0 < b_I.size(0); idx0++) {
			for (int idx1{ 0 }; idx1 < b_I.size(1); idx1++) {
				b_I[idx0 + b_I.size(0) * idx1] = threshold_image.at<uchar>(idx0, idx1);
			}
		}
		imedge_2d(b_I, 0, binary_I);

		for (int idx0{ 0 }; idx0 < b_I.size(0); idx0++) {
			for (int idx1{ 0 }; idx1 < b_I.size(1); idx1++)
			{
				if (binary_I[idx0 + b_I.size(0) * idx1])
				{
					threshold_image.at<uchar>(idx0, idx1) = 255;
				}
				else
				{
					threshold_image.at<uchar>(idx0, idx1) = 0;
				}
			}
		}
		findContours(threshold_image, contours, hierarchy, cv::RETR_LIST, cv::CHAIN_APPROX_NONE, cv::Point(0, 0));
	}
	if (image_process_method == CANNY_Method)
	{
		threshold_image = ori_image_mat.clone();
		if (threshold_image.channels() == 3)
		{
			cv::cvtColor(threshold_image, threshold_image, CV_BGR2GRAY);
		}
		coder::array<unsigned char, 2U> binary_I;
		binary_I.set_size(threshold_image.rows, threshold_image.cols);
		coder::array<boolean_T, 2U> b_I;
		b_I.set_size(threshold_image.rows, threshold_image.cols);
		for (int idx0{ 0 }; idx0 < b_I.size(0); idx0++) {
			for (int idx1{ 0 }; idx1 < b_I.size(1); idx1++) {
				b_I[idx0 + b_I.size(0) * idx1] = threshold_image.at<uchar>(idx0, idx1);
			}
		}
		imedge_2d(b_I, 1, binary_I);

		for (int idx0{ 0 }; idx0 < b_I.size(0); idx0++) {
			for (int idx1{ 0 }; idx1 < b_I.size(1); idx1++) 
			{
				if (binary_I[idx0 + b_I.size(0) * idx1])
				{
					threshold_image.at<uchar>(idx0, idx1) = 255;
				}
				else
				{
					threshold_image.at<uchar>(idx0, idx1) = 0;
				}
			}
		}
		findContours(threshold_image, contours, hierarchy, cv::RETR_LIST, cv::CHAIN_APPROX_NONE, cv::Point(0, 0));
	}
	if (image_process_method == Prewitt_Method)
	{
		threshold_image = ori_image_mat.clone();
		if (threshold_image.channels() == 3)
		{
			cv::cvtColor(threshold_image, threshold_image, CV_BGR2GRAY);
		}
		coder::array<unsigned char, 2U> binary_I;
		binary_I.set_size(threshold_image.rows, threshold_image.cols);
		coder::array<boolean_T, 2U> b_I;
		b_I.set_size(threshold_image.rows, threshold_image.cols);
		for (int idx0{ 0 }; idx0 < b_I.size(0); idx0++) {
			for (int idx1{ 0 }; idx1 < b_I.size(1); idx1++) {
				b_I[idx0 + b_I.size(0) * idx1] = threshold_image.at<uchar>(idx0, idx1);
			}
		}
		imedge_2d(b_I, 2, binary_I);

		for (int idx0{ 0 }; idx0 < b_I.size(0); idx0++) {
			for (int idx1{ 0 }; idx1 < b_I.size(1); idx1++)
			{
				if (binary_I[idx0 + b_I.size(0) * idx1])
				{
					threshold_image.at<uchar>(idx0, idx1) = 255;
				}
				else
				{
					threshold_image.at<uchar>(idx0, idx1) = 0;
				}
			}
		}	
		findContours(threshold_image, contours, hierarchy, cv::RETR_LIST, cv::CHAIN_APPROX_NONE, cv::Point(0, 0));
	}
	if (image_process_method == Roberts_Method)
	{
		threshold_image = ori_image_mat.clone();
		if (threshold_image.channels() == 3)
		{
			cv::cvtColor(threshold_image, threshold_image, CV_BGR2GRAY);
		}
		coder::array<unsigned char, 2U> binary_I;
		binary_I.set_size(threshold_image.rows, threshold_image.cols);
		coder::array<boolean_T, 2U> b_I;
		b_I.set_size(threshold_image.rows, threshold_image.cols);
		for (int idx0{ 0 }; idx0 < b_I.size(0); idx0++) {
			for (int idx1{ 0 }; idx1 < b_I.size(1); idx1++) {
				b_I[idx0 + b_I.size(0) * idx1] = threshold_image.at<uchar>(idx0, idx1);
			}
		}
		imedge_2d(b_I, 3, binary_I);

		for (int idx0{ 0 }; idx0 < b_I.size(0); idx0++) {
			for (int idx1{ 0 }; idx1 < b_I.size(1); idx1++)
			{
				if (binary_I[idx0 + b_I.size(0) * idx1])
				{
					threshold_image.at<uchar>(idx0, idx1) = 255;
				}
				else
				{
					threshold_image.at<uchar>(idx0, idx1) = 0;
				}
			}
		}
		findContours(threshold_image, contours, hierarchy, cv::RETR_LIST, cv::CHAIN_APPROX_NONE, cv::Point(0, 0));
	}
	else if (image_process_method == OTSU_Method)
	{
		threshold(ori_image_mat, threshold_image, 0, 255.0, cv::THRESH_OTSU | cv::THRESH_BINARY_INV);

		findContours(threshold_image, contours, hierarchy, cv::RETR_LIST, cv::CHAIN_APPROX_NONE, cv::Point(0, 0));
	}
	else if (image_process_method == ADAPTIVE_THRESH_Method)
	{
		adaptiveThreshold(ori_image_mat, threshold_image, 255, cv::ADAPTIVE_THRESH_MEAN_C, cv::THRESH_BINARY, 5, 0);
		findContours(threshold_image, contours, hierarchy, cv::RETR_LIST, cv::CHAIN_APPROX_NONE, cv::Point(0, 0));
	}
	else if (image_process_method == Block_OTSU)
	{
		threshold(ori_image_mat, threshold_image, 0, 255.0, cv::THRESH_OTSU | cv::THRESH_BINARY_INV);

		findContours(threshold_image, contours, hierarchy, cv::RETR_LIST, cv::CHAIN_APPROX_NONE, cv::Point(0, 0));
	}
	//cv::namedWindow("canny", cv::WINDOW_NORMAL);
	//cv::imshow("canny", threshold_image);
	//cv::waitKey(0);
	return true;
}

float ImageDetectMethod::Contours_arc(std::vector<cv::Point> C , cv::Point2f center)
{
	if (C.size() < 2)
	{
		return 0;
	}
	std::vector<float> arc_vec;
	for (int ii = 0; ii < C.size(); ii++)
	{
		arc_vec.push_back(atan2f((float)C[ii].x - center.x, (float)C[ii].y - center.y));
	}
	std::sort(arc_vec.begin(), arc_vec.end());
	double arc_max = abs(arc_vec[0] - arc_vec[arc_vec.size() - 1]);
	arc_max = arc_max < (6.283185307179586 - arc_max) ? arc_max : (6.283185307179586 - arc_max);
	for (int ii = 0; ii < C.size() - 1; ii++)
	{
		double L_now = abs(arc_vec[ii + 1] - arc_vec[ii]);
		L_now = L_now < (6.283185307179586 - L_now) ? L_now : (6.283185307179586 - L_now);
		arc_max = arc_max < L_now ? L_now : arc_max;
	}
	return (6.283185307179586 - arc_max) / 6.283185307179586 * 360.0;
}
float ImageDetectMethod::Contours_Length(std::vector<cv::Point> C)
{
	if (C.size() < 2)
	{
		return 0;
	}
	double L_max = sqrt(pow(C[0].x - C[C.size() - 1].x, 2) + pow(C[0].y - C[C.size() - 1].y, 2));
	double L_sum = L_max;
	for (int ii = 0; ii < C.size() - 1; ii++)
	{
		double L_now = sqrt(pow(C[ii].x - C[ii + 1].x, 2) + pow(C[ii].y - C[ii + 1].y, 2));
		L_max = L_max < L_now ? L_now : L_max;
		L_sum += L_now;
	}
	return L_sum - L_max;
}
bool ImageDetectMethod::FilterEllipseContours(const std::vector<std::vector<cv::Point>>& contours
	, int min_radius_pixel, int max_radius_pixel,
	float ellipse_error_pixel, 
	QList<QList<float>>& ellipse_pars, double max_aspect_ratio, int min_points,int min_contour_num,
	int min_arc)
{
	//ellipse_error_pixel = 0.5;

	cv::Mat ori_image = 255 * cv::Mat::zeros(3648, 5472, CV_8U);
	min_points = min_points < 6 ? 6 : min_points;
	ellipse_pars.clear();
	for (int i = 0; i < contours.size(); i++)
	{
		//椭圆最小6点拟合
		int count = contours[i].size();
		if (count < min_points)
		{
			continue;
		}

		//周长,面积
		double length = arcLength(contours[i], true);
		double area = contourArea(contours[i]);
		double perimeter = cv::arcLength(contours[i], true);
		std::vector<cv::Point> contours_ploy;
		cv::approxPolyDP(contours[i], contours_ploy, 0.02 * perimeter, true);
		int CornerNum = contours_ploy.size();
		if (CornerNum < min_contour_num)
		{
			continue;
		}
		cv::RotatedRect rRect = fitEllipseAMS(contours[i]);
		if (isnan(rRect.center.x) || isinf(rRect.center.x))
		{
			continue;
		}
		if (isnan(rRect.center.y) || isinf(rRect.center.y))
		{
			continue;
		}
		if (isnan(rRect.size.width) || isinf(rRect.size.width))
		{
			continue;
		}
		if (isnan(rRect.size.height) || isinf(rRect.size.height))
		{
			continue;
		}

		//半径大小判断
		if (2 * min_radius_pixel >= rRect.size.width || 2 * max_radius_pixel < +rRect.size.height)
		{
			continue;
		}

		////////////////长短轴比例判断
		if (rRect.size.height / rRect.size.width > max_aspect_ratio)
		{
			continue;
		}

		////////////////最小二乘拟合误差判断，半径差值


		float whole_error = 0;
		for (int j = 0; j < count; j++)
		{
			whole_error += ErrorDROfEllipseFit(rRect.center.x, rRect.center.y, rRect.size.width * 0.5, rRect.size.height * 0.5, rRect.angle * M_PI / 180,
				float(contours[i][j].x), float(contours[i][j].y));
		}
		float aver_error = whole_error / count;

		if (aver_error > ellipse_error_pixel * 1.5)
		{
			continue;
		}
		if (Contours_arc(contours[i], rRect.center)< min_arc)
		{
			continue;
		}
		////二次优化
		//double* coff_now = new double[5];
		//coff_now[0] = rRect.center.x;
		//coff_now[1] = rRect.center.y;
		//coff_now[2] = rRect.size.width / 2.0;
		//coff_now[3] = rRect.size.height / 2.0;
		//coff_now[4] = rRect.angle * M_PI / 180;
		//ceres::Problem problem;
		//for (int jj = 0; jj < contours[i].size(); jj++)
		//{
		//	Eigen::Vector2d ellipse_point(contours[i][jj].x, contours[i][jj].y);
		//	ceres::CostFunction* cost_function = Fit_ellipse_ceres::Create(ellipse_point);
		//	ceres::LossFunction* loss_function = new ceres::CauchyLoss(1);
		//	problem.AddResidualBlock(cost_function, loss_function, coff_now);
		//}
		//ceres::Solver::Summary summary;
		//ceres::Solver::Options options;
		//options.minimizer_progress_to_stdout = false;
		//options.max_num_iterations = 1000;
		//options.function_tolerance = 1e-9;
		//options.gradient_tolerance = 1e-9;
		//options.parameter_tolerance = 1e-9;
		//ceres::Solve(options, &problem, &summary);


		//cv::RotatedRect rRect_ceres;
		//rRect_ceres.center.x = coff_now[0];
		//rRect_ceres.center.y = coff_now[1];
		//rRect_ceres.angle = coff_now[4] / M_PI * 180;
		//rRect_ceres.size.width = 2 * coff_now[2];
		//rRect_ceres.size.height = 2 * coff_now[3];
		//float whole_error_ceres = 0;
		//for (int j = 0; j < count; j++)
		//{
		//	whole_error_ceres += ErrorDROfEllipseFit(rRect_ceres.center.x, rRect_ceres.center.y, rRect_ceres.size.width * 0.5, rRect_ceres.size.height * 0.5, rRect_ceres.angle * M_PI / 180,
		//		float(contours[i][j].x), float(contours[i][j].y));
		//}
		//float aver_error_ceres = whole_error_ceres / count;

		////if (aver_error_ceres > ellipse_error_pixel)
		////{
		////	continue;
		////}
		//QList<float> ellipse_par_list;
		//ellipse_par_list << rRect_ceres.center.x;
		//ellipse_par_list << rRect_ceres.center.y;
		//ellipse_par_list << rRect_ceres.size.width * 0.5;
		//ellipse_par_list << rRect_ceres.size.height * 0.5;
		//ellipse_par_list << rRect_ceres.angle;
		//ellipse_par_list << aver_error_ceres;
		//ellipse_par_list << i;

		QList<float> ellipse_par_list;
		ellipse_par_list << rRect.center.x;
		ellipse_par_list << rRect.center.y;
		ellipse_par_list << rRect.size.width * 0.5;
		ellipse_par_list << rRect.size.height * 0.5;
		ellipse_par_list << rRect.angle * M_PI / 180;
		ellipse_par_list << aver_error;
		ellipse_par_list << i;
		ellipse_pars.append(ellipse_par_list);
	}
	return true;
}

float ImageDetectMethod::ErrorDROfEllipseFit(float center_x, float center_y, float ellipse_a, float ellipse_b,
	float ellipse_angle_in_pi, float x, float y)
{
	//坐标转换，图像坐标系转椭圆坐标系      x^2/a^2 + y^2/b^2 =1
	float tr_x = (x - center_x) * cos(ellipse_angle_in_pi) + (y - center_y) * sin(ellipse_angle_in_pi);
	float tr_y = -(x - center_x) * sin(ellipse_angle_in_pi) + (y - center_y) * cos(ellipse_angle_in_pi);

	//计算拟合误差，半径差作为拟合误差
	float alfa = atan2(tr_y, tr_x);

	float r = ellipse_a * ellipse_b / sqrt(ellipse_a * ellipse_a * sin(alfa) * sin(alfa) + ellipse_b * ellipse_b * cos(alfa) * cos(alfa));
	float delta_r = sqrt(tr_x * tr_x + tr_y * tr_y) - r;

	return abs(delta_r);
}

float ImageDetectMethod::LeastSquareErrorOfEllipseFit(float center_x, float center_y, float ellipse_a, float ellipse_b,
	float ellipse_angle_in_pi, float x, float y)
{
	//坐标转换，图像坐标系转椭圆坐标系          x^2/a^2 + y^2/b^2 =1
	float tr_x = (x - center_x) * cos(ellipse_angle_in_pi) + (y - center_y) * sin(ellipse_angle_in_pi);
	float tr_y = -(x - center_x) * sin(ellipse_angle_in_pi) + (y - center_y) * cos(ellipse_angle_in_pi);

	//计算拟合误差，最小二乘拟合误差
	float e = tr_x * tr_x / ellipse_a / ellipse_a + tr_y * tr_y / ellipse_b / ellipse_b - 1;

	return e * e;
}

bool ImageDetectMethod::FilterEllipseContoursForCodePoint(const cv::Mat& image_mat, float ratio_k, float ratio_k1, float ratio_k2,
	QList<QList<float>>& ellipse_pars,
	float delta_Mt, float fore_stdDev, float back_stdDev)//code uncode的区分，C
{
	//进一步筛选出用编码点解码的点
	//灰度准则
	for (int i = 0; i < ellipse_pars.size(); i++)
	{
		float out_foreground_mean = 0;
		float out_background_mean = 0;
		float out_foreground_std = 0;
		float out_background_std = 0;
		bool is_gray_judge = EllipseGrayJudgeForCodePoint(image_mat, ellipse_pars[i][0], ellipse_pars[i][1],
			ellipse_pars[i][2], ellipse_pars[i][3], ellipse_pars[i][4], ratio_k, out_foreground_mean, out_background_mean
			, out_foreground_std, out_background_std, delta_Mt, fore_stdDev, back_stdDev);

		if (is_gray_judge == false)
		{
			ellipse_pars.removeAt(i);
			i--;
		}
	}
	//剔除小的编码带影响
	///////筛选出标志点进行解码,位置关系，剔除外圆环和小编码带
	for (int i = 0; i < ellipse_pars.size() - 1; i++)
	{
		for (int j = i + 1; j < ellipse_pars.size(); j++)
		{
			float x1 = ellipse_pars[i][0];
			float y1 = ellipse_pars[i][1];
			float x2 = ellipse_pars[j][0];
			float y2 = ellipse_pars[j][1];
			float length_of_2_points = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));

			if (length_of_2_points < std::min(ellipse_pars[i][3], ellipse_pars[j][3]))
			{
				if (ellipse_pars[i][3] > ellipse_pars[j][3])
				{
					ellipse_pars.removeAt(i);
					i--;
					break;
				}
				else
				{
					ellipse_pars.removeAt(j);
					j--;
				}
			}
			else
				if (length_of_2_points < std::min(ellipse_pars[i][3], ellipse_pars[j][3]) * ratio_k2)
				{
					if (ellipse_pars[i][5] / sqrt(ellipse_pars[i][2] * ellipse_pars[i][3]) > ellipse_pars[j][5] / sqrt(ellipse_pars[j][2] * ellipse_pars[j][3]))
					{
						ellipse_pars.removeAt(i);
						i--;
						break;
					}
					else
					{
						ellipse_pars.removeAt(j);
						j--;
					}
				}
		}
	}
	return true;
}
bool ImageDetectMethod::EllipseGrayJudgeForPointCSI_is2Circle(
	QList<QList<float>> ellipse_pars_all, QList<float> ellipse_pars_now, float rati_k)
{
	bool is_exist = false;
	for (int i = 0; i < ellipse_pars_all.size() - 1; i++)
	{
		if (ellipse_pars_all[i][6] == ellipse_pars_now[6])
		{
			continue;
		}
		float x1 = ellipse_pars_all[i][0];
		float y1 = ellipse_pars_all[i][1];
		float x2 = ellipse_pars_now[0];
		float y2 = ellipse_pars_now[1];
		float length_of_2_points = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));

		if (length_of_2_points < std::min(ellipse_pars_all[i][3], ellipse_pars_now[3]))
		{
			double r1 = ellipse_pars_all[i][3] < ellipse_pars_now[3] ? ellipse_pars_all[i][3] : ellipse_pars_now[3];
			double r2 = ellipse_pars_all[i][3] > ellipse_pars_now[3] ? ellipse_pars_all[i][3] : ellipse_pars_now[3];
			double r_now = r1 / r2;
			if (r_now > rati_k * 0.5 && r_now < rati_k * 1.5)
			{
				is_exist = true;
				break;
			}
		}
	}
	return is_exist;
}
bool ImageDetectMethod::EllipseGrayJudgeForPointCSI(const cv::Mat& image_mat, float center_x, float center_y,
	float ellipse_a, float ellipse_b, float angle_in_pi, float ratio_k,
	float& out_std)
{
	float ellipse_a2;
	float ellipse_b2;

	if (ratio_k > 1)
	{
		ellipse_a2 = ellipse_a * ratio_k;
		ellipse_b2 = ellipse_b * ratio_k;
	}
	else
	{
		ellipse_a2 = ellipse_a;
		ellipse_b2 = ellipse_b;
		ellipse_a = ellipse_a2 * ratio_k;
		ellipse_b = ellipse_b2 * ratio_k;
	}

	int i_top = 0;
	int i_bottom = 0;
	int j_left = 0;
	int j_right = 0;
	if (int(center_y - ellipse_b2) < 0)
	{
		i_top = 0;
	}
	else i_top = int(center_y - ellipse_b2);
	if (int(center_y + ellipse_b2) > image_mat.rows - 1)
	{
		i_bottom = image_mat.rows - 1;
	}
	else i_bottom = int(center_y + ellipse_b2);
	if (int(center_x - ellipse_b2) < 0)
	{
		j_left = 0;
	}
	else j_left = int(center_x - ellipse_b2);
	if (int(center_x + ellipse_b2) > image_mat.cols - 1)
	{
		j_right = image_mat.cols - 1;
	}
	else j_right = int(center_x + ellipse_b2);

	float foreground_median = 0;
	float background_median = 0;

	QList<float> background_value_list;
	QList<float> foreground_value_list;
	QList<float> all_value_list;
	QList<float> whole_value_list;
	cv::Mat cropped_image;
	cv::Mat PROCE_img;
	for (int i = i_top; i <= i_bottom; i++)
	{
		for (int j = j_left; j <= j_right; j++)
		{
			const uchar* _image_mat_ptr = image_mat.ptr<uchar>(i);


			float tr_x = (j - center_x) * cos(angle_in_pi) + (i - center_y) * sin(angle_in_pi);
			float tr_y = -(j - center_x) * sin(angle_in_pi) + (i - center_y) * cos(angle_in_pi);

			if (tr_x * tr_x / int(ellipse_a2) / int(ellipse_a2) + tr_y * tr_y / int(ellipse_b2) / int(ellipse_b2) < 1 &&
				tr_x * tr_x / ellipse_a / ellipse_a + tr_y * tr_y / ellipse_b / ellipse_b > 1)
			{
				background_value_list.append(float(_image_mat_ptr[j]));
				//std::cout << 0 << "\t" << i << "\t" << j << "\t" << float(_image_mat_ptr[j]) << "\n";
			}
			else if (tr_x * tr_x / int(ellipse_a) / int(ellipse_a) + tr_y * tr_y / int(ellipse_b) / int(ellipse_b) < 1)
			{
				foreground_value_list.append(float(_image_mat_ptr[j]));
				//std::cout << 1 << "\t" << i << "\t" << j << "\t" << float(_image_mat_ptr[j]) << "\n";
			}
		}
	}
	int inter_1 = 1, inter_2 = 1;
	if (ratio_k * ratio_k - 1 > 1)
	{
		inter_1 = round(ratio_k * ratio_k - 1);
		inter_1 = inter_1 < 1 ? 1 : inter_1;
		inter_2 = 1;
	}
	else
	{
		inter_2 = round(1.0 / (ratio_k * ratio_k - 1));
		inter_2 = inter_2 < 1 ? 1 : inter_2;
		inter_1 = 1;
	}
	double mean_value = 0;
	double mean_value_number = 0;
	for (int ii = 0; ii < background_value_list.size(); ii = ii + inter_1)
	{
		all_value_list.push_back(background_value_list[ii]);
		mean_value += background_value_list[ii];
		mean_value_number += 1.0;
	}
	for (int ii = 0; ii < foreground_value_list.size(); ii = ii = ii + inter_2)
	{
		all_value_list.push_back(foreground_value_list[ii]);
		mean_value += foreground_value_list[ii];
		mean_value_number += 1.0;
	}
	mean_value /= mean_value_number;

	double std_value = 0;
	for (int ii = 0; ii < all_value_list.size(); ii++)
	{
		std_value += pow(all_value_list[ii] - mean_value, 2);
	}

	std_value /= mean_value_number;
	std_value = sqrt(std_value);
	out_std = std_value;

	return true;
}
bool ImageDetectMethod::EllipseGrayJudgeForCodePoint(const cv::Mat& image_mat, float center_x, float center_y,
	float ellipse_a, float ellipse_b, float angle_in_pi, float ratio_k,
	float& out_foreground_mean, float& out_background_mean, float& out_foreground_stdDev, float& out_background_stdDev,
	float delta_Mt, float fore_stdDev, float back_stdDev)
{
	float ellipse_a2 = ellipse_a * ratio_k;
	float ellipse_b2 = ellipse_b * ratio_k;

	int i_top = 0;
	int i_bottom = 0;
	int j_left = 0;
	int j_right = 0;
	if (int(center_y - ellipse_b2) < 0)
	{
		i_top = 0;
	}
	else i_top = int(center_y - ellipse_b2);
	if (int(center_y + ellipse_b2) > image_mat.rows - 1)
	{
		i_bottom = image_mat.rows - 1;
	}
	else i_bottom = int(center_y + ellipse_b2);
	if (int(center_x - ellipse_b2) < 0)
	{
		j_left = 0;
	}
	else j_left = int(center_x - ellipse_b2);
	if (int(center_x + ellipse_b2) > image_mat.cols - 1)
	{
		j_right = image_mat.cols - 1;
	}
	else j_right = int(center_x + ellipse_b2);

	float foreground_mean = 0;
	float background_mean = 0;
	float whole_mean = 0;
	float foreground_stdDev = 0;
	float background_stdDev = 0;
	float whole_stdDev = 0;
	int forground_num = 0;
	int background_num = 0;
	int whole_num = 0;

	QList<float> background_value_list;
	QList<float> foreground_value_list;
	QList<float> whole_value_list;
	cv::Mat cropped_image;
	cv::Mat PROCE_img;

	for (int i = i_top; i <= i_bottom; i++)
	{
		for (int j = j_left; j <= j_right; j++)
		{
			const uchar* _image_mat_ptr = image_mat.ptr<uchar>(i);
			//坐标转换图像坐标系转椭圆坐标系


			float tr_x = (j - center_x) * cos(angle_in_pi) + (i - center_y) * sin(angle_in_pi);
			float tr_y = -(j - center_x) * sin(angle_in_pi) + (i - center_y) * cos(angle_in_pi);
			//float tr_x = (j - center_x) * sin(angle_in_pi) - (i - center_y) * cos(angle_in_pi);
			//float tr_y = (j - center_x) * cos(angle_in_pi) + (i - center_y) * sin(angle_in_pi);

			if (tr_x * tr_x / int(ellipse_a2) / int(ellipse_a2) + tr_y * tr_y / int(ellipse_b2) / int(ellipse_b2) < 1 &&
				tr_x * tr_x / ellipse_a / ellipse_a + tr_y * tr_y / ellipse_b / ellipse_b > 1)
			{
				background_mean += float(_image_mat_ptr[j]);
				background_value_list.append(float(_image_mat_ptr[j]));
				background_num++;

				whole_mean += float(_image_mat_ptr[j]);
				whole_value_list.append(float(_image_mat_ptr[j]));
				whole_num++;
			}
			else if (tr_x * tr_x / int(ellipse_a) / int(ellipse_a) + tr_y * tr_y / int(ellipse_b) / int(ellipse_b) < 1)
			{
				foreground_mean += float(_image_mat_ptr[j]);
				foreground_value_list.append(float(_image_mat_ptr[j]));
				forground_num++;

				whole_mean += float(_image_mat_ptr[j]);
				whole_value_list.append(float(_image_mat_ptr[j]));
				whole_num++;
			}
		}
	}
	foreground_mean = foreground_mean / forground_num;
	background_mean = background_mean / background_num;
	whole_mean = whole_mean / whole_num;
	for (int i = 0; i < background_value_list.size(); i++)
	{
		background_stdDev += (background_value_list[i] - background_mean) * (background_value_list[i] - background_mean);
	}
	for (int i = 0; i < foreground_value_list.size(); i++)
	{
		foreground_stdDev += (foreground_value_list[i] - foreground_mean) * (foreground_value_list[i] - foreground_mean);
	}
	for (int i = 0; i < whole_value_list.size(); i++)
	{
		whole_stdDev += (whole_value_list[i] - whole_mean) * (whole_value_list[i] - whole_mean);
	}
	foreground_stdDev = sqrt(foreground_stdDev / forground_num);
	background_stdDev = sqrt(background_stdDev / background_num);
	whole_stdDev = sqrt(whole_stdDev / whole_num);
	out_foreground_mean = foreground_mean;
	out_background_mean = background_mean;
	out_foreground_stdDev = foreground_stdDev;
	out_background_stdDev = background_stdDev;
	if (abs(foreground_mean - background_mean) < delta_Mt)
	{
		return false;
	}
	if (foreground_stdDev > fore_stdDev)
	{
		return false;
	}
	if (background_stdDev > back_stdDev)
	{
		return false;
	}

	return true;
}

int* ImageDetectMethod::ReturnDefualtIdArray(int& array_size, CodePointBitesType code_bites_type /*=CodeBites15*/)
{
	switch (code_bites_type)
	{
	case CodeBites15:
		array_size = 429;
		//array_size = sizeof(static_15_8_code_id_false) / sizeof(static_15_8_code_id_false[0]);

		return &static_15_8_code_id_true[0];
		break;
	}

	return NULL;
}

bool ImageDetectMethod::Decoding20140210(const cv::Mat& image_mat, int& out_put_code_id, float center_x, float center_y,
	float ellipse_a, float ellipse_b, float angle_in_pi,
	float ratio_k1 /*=2.4*/, float ratio_k2 /*=4 */,
	MarkPointColorType color_type /*= BlackDownWhiteUp*/, CodePointBitesType code_bites_type /*=CodeBites15*/, double thread_value/* = 0*/)
{
	// 	if (center_x>679&&center_x<721&&center_y>568&&center_y<590)
	// 	{
	// 		int a =0;
	// 	}
	float ellipse_a2 = ellipse_a * ((ratio_k2 + ratio_k1) / 2.0 - (ratio_k2 - ratio_k1) / 2.0 * 0.75);
	float ellipse_b2 = ellipse_b * ((ratio_k2 + ratio_k1) / 2.0 - (ratio_k2 - ratio_k1) / 2.0 * 0.75);
	float ellipse_a4 = ellipse_a * ((ratio_k2 + ratio_k1) / 2.0 + (ratio_k2 - ratio_k1) / 2.0 * 0.75);
	float ellipse_b4 = ellipse_b * ((ratio_k2 + ratio_k1) / 2.0 + (ratio_k2 - ratio_k1) / 2.0 * 0.75);

	int i_top = 0;
	int i_bottom = 0;
	int j_left = 0;
	int j_right = 0;
	if (int(center_y - ellipse_b4) < 0)
	{
		return false;
	}
	else i_top = int(center_y - ellipse_b4);
	if (int(center_y + ellipse_b4) > image_mat.rows - 1)
	{
		return false;
	}
	else i_bottom = int(center_y + ellipse_b4);
	if (int(center_x - ellipse_b4) < 0)
	{
		return false;
	}
	else j_left = int(center_x - ellipse_b4);
	if (int(center_x + ellipse_b4) > image_mat.cols - 1)
	{
		return false;
	}
	else j_right = int(center_x + ellipse_b4);

	int coed_bites_num;
	switch (code_bites_type)
	{
	case CodeBites15:
		coed_bites_num = 15;
		break;
	case CodeBites12:
		coed_bites_num = 12;
		break;
	}

	QList<QPointF> ellipse_ring_points;
	QList<QPointF> ellipse_ring_points_for_cal_1;
	QList<QPointF> ellipse_ring_points_for_cal_2;
	QList<int> ellipse_ring_gray_value;

	//将圆环分成15段，每段取20个点的值
	int per_num_point = 30;

	for (int i = 0; i < coed_bites_num * per_num_point; i++)
	{
		float x = cos(float(i) / coed_bites_num / per_num_point * 2 * M_PI);
		float y = sin(float(i) / coed_bites_num / per_num_point * 2 * M_PI);

		float ell_cor_r2_x = ellipse_a2 * x;
		float ell_cor_r2_y = ellipse_b2 * y;
		float ell_cor_r3_x = ellipse_a4 * x;
		float ell_cor_r3_y = ellipse_b4 * y;


		//float tr_x = (j - center_x) * sin(angle_in_pi) - (i - center_y) * cos(angle_in_pi);
		//float tr_y = (j - center_x) * cos(angle_in_pi) + (i - center_y) * sin(angle_in_pi);

		float car_cor_r2_x = ell_cor_r2_x * cos(angle_in_pi) - ell_cor_r2_y * sin(angle_in_pi) + center_x;
		float car_cor_r2_y = ell_cor_r2_x * sin(angle_in_pi) + ell_cor_r2_y * cos(angle_in_pi) + center_y;
		float car_cor_r3_x = ell_cor_r3_x * cos(angle_in_pi) - ell_cor_r3_y * sin(angle_in_pi) + center_x;
		float car_cor_r3_y = ell_cor_r3_x * sin(angle_in_pi) + ell_cor_r3_y * cos(angle_in_pi) + center_y;
		//float car_cor_r2_x = ell_cor_r2_x * sin(angle_in_pi) + ell_cor_r2_y * cos(angle_in_pi) + center_x;
		//float car_cor_r2_y = -ell_cor_r2_x * cos(angle_in_pi) + ell_cor_r2_y * sin(angle_in_pi) + center_y;
		//float car_cor_r3_x = ell_cor_r3_x * sin(angle_in_pi) + ell_cor_r3_y * cos(angle_in_pi) + center_x;
		//float car_cor_r3_y = -ell_cor_r3_x * cos(angle_in_pi) + ell_cor_r3_y * sin(angle_in_pi) + center_y;
		QPoint point1 = QPoint(int(car_cor_r2_x), int(car_cor_r2_y));
		QPoint point2 = QPoint(int(car_cor_r3_x), int(car_cor_r3_y));

		QList<int> gray_list = GetALineGrayList(image_mat, point1, point2);
		//std::cout << point1.x() << "\t" << point1.y() << "\t" << point2.x() << "\t" << point2.y() << "\n";
		int mid_value = MIdValue(gray_list);
		ellipse_ring_points_for_cal_1.append(QPointF(car_cor_r2_x, car_cor_r2_y));
		ellipse_ring_points_for_cal_2.append(QPointF(car_cor_r3_x, car_cor_r3_y));
		ellipse_ring_points.append(QPointF(x, y));
		ellipse_ring_gray_value.append(mid_value);
	}
	QList<int> ellipse_ring_gray_value_new;
	for (int ii = 2; ii < ellipse_ring_gray_value.size() - 2; ii++)
	{
		QList<int> mid_list = { ellipse_ring_gray_value[ii - 2],ellipse_ring_gray_value[ii - 1], ellipse_ring_gray_value[ii], ellipse_ring_gray_value[ii + 2],
		ellipse_ring_gray_value[ii + 2] };
		int mid_value_list = MIdValue(mid_list);
		if (ii == 2)
		{
			ellipse_ring_gray_value_new.append(mid_value_list);
			ellipse_ring_gray_value_new.append(mid_value_list);
			ellipse_ring_gray_value_new.append(mid_value_list);
		}
		else if (ii == ellipse_ring_gray_value.size() - 3)
		{
			ellipse_ring_gray_value_new.append(mid_value_list);
			ellipse_ring_gray_value_new.append(mid_value_list);
			ellipse_ring_gray_value_new.append(mid_value_list);
		}
		else
		{
			ellipse_ring_gray_value_new.append(mid_value_list);
		}
	}
	ellipse_ring_gray_value = ellipse_ring_gray_value_new;
	////找单位圆编码的边界，一维卷积模板半宽为3
	int half_length = per_num_point;
	int edge_index = 0;
	float ans_max = 0;
	for (int i = 0; i < ellipse_ring_points.size(); i++)
	{
		float ans = 0;
		for (int j = 0; j < half_length; j++)
		{
			int index1 = i + j;
			if (index1 >= ellipse_ring_gray_value.size())
			{
				index1 = index1 % ellipse_ring_gray_value.size();
			}
			int index2 = i + j + half_length;
			if (index2 >= ellipse_ring_gray_value.size())
			{
				index2 = index2 % ellipse_ring_gray_value.size();
			}
			ans += -ellipse_ring_gray_value.at(index1) + ellipse_ring_gray_value.at(index2);
		}
		if (abs(ans) > ans_max)
		{
			ans_max = abs(ans);
			edge_index = i + half_length;
		}
	}

	if (ans_max / half_length < 20)
	{
		return false;
	}
	//////////解编码，24度顺时针

	QList<int> code_in_2;
	double threshold_value = 0 /*= 150*/;       //阈值,要修改
	for (int i = edge_index - half_length; i < edge_index + half_length; i++)
	{
		int index1 = i;
		if (index1 >= ellipse_ring_gray_value.size())
		{
			index1 = index1 % ellipse_ring_gray_value.size();
		}
		threshold_value += ellipse_ring_gray_value[index1];
	}
	threshold_value /= half_length * 2;
	double bit_0_value = 0;
	double bit_1_value = 0;
	int bit_0_number = 0;
	int bit_1_number = 0;
	for (int i = 0; i < coed_bites_num; i++)
	{
		double mean_gray = 0;
		int num = 0;
		for (int j = 0; j < per_num_point; j++)
		{
			int index = edge_index + i * per_num_point + j + 1;
			if (index >= ellipse_ring_gray_value.size())
			{
				index = index % ellipse_ring_gray_value.size();
			}
			mean_gray += ellipse_ring_gray_value[index];
		}
		mean_gray /= per_num_point;

		if (color_type == BlackDownWhiteUp)
		{
			if (mean_gray > threshold_value)
			{
				bit_1_value += mean_gray;
				bit_1_number++;
				code_in_2.append(1);
			}
			else
			{
				bit_0_value += mean_gray;
				bit_0_number++;
				code_in_2.append(0);
			}
		}
		else
		{
			if (mean_gray > threshold_value)
			{
				bit_0_value += mean_gray;
				bit_0_number++;
				code_in_2.append(0);
			}
			else
			{
				bit_1_value += mean_gray;
				bit_1_number++;
				code_in_2.append(1);
			}
		}
	}
	bit_1_value /= (double)bit_1_number;
	bit_0_value /= (double)bit_0_number;
	if (abs(bit_1_value- bit_0_value) < thread_value)
	{
		return false;
	}
	int code_id;
	QList<int> code_out_2;
	CalculateRealCodeID20140210(code_in_2, code_out_2, code_id);

	out_put_code_id = code_id;

	return true;
}

bool ImageDetectMethod::CalculateRealCodeID20140210(QList<int> in_put_code_list, QList<int>& out_put_code_list, int& out_put_code_ID)
{
	out_put_code_ID = Change2To10(in_put_code_list);
	out_put_code_list = in_put_code_list;

	int n = in_put_code_list.size();

	for (int i = 1; i <= n - 1; i++)
	{
		QList<int> new_code_list;
		int new_id;
		for (int j = 0; j < n; j++)
		{
			if (i + j <= n - 1)
			{
				new_code_list.append(in_put_code_list.at(i + j));
			}
			else
			{
				new_code_list.append(in_put_code_list.at(i + j - n));
			}
		}
		new_id = Change2To10(new_code_list);
		if (out_put_code_ID > new_id)
		{
			out_put_code_ID = new_id;
			out_put_code_list = new_code_list;
		}
	}
	return true;
}

bool ImageDetectMethod::UncodePointCheck(const cv::Mat& image_mat, float center_x, float center_y, float ellipse_a, float ellipse_b,
	float angle_in_pi, float ratio_k /*=2*/,
	MarkPointColorType color_type /*= BlackDownWhiteUp*/, CodePointBitesType code_bites_type /*=CodeBites15*/)//非编码点的判断，c
{
	double ellipse_a1 = ellipse_a * 1.1;
	double ellipse_b1 = ellipse_b * 1.1;
	double ellipse_a2 = ellipse_a * (ratio_k - 0.1);
	double ellipse_b2 = ellipse_b * (ratio_k - 0.1);

	int i_top = 0;
	int i_bottom = 0;
	int j_left = 0;
	int j_right = 0;
	if (int(center_y - ellipse_b2) < 0)
	{
		return false;
	}
	else i_top = int(center_y - ellipse_b2);
	if (int(center_y + ellipse_b2) > image_mat.rows - 1)
	{
		return false;
	}
	else i_bottom = int(center_y + ellipse_b2);
	if (int(center_x - ellipse_b2) < 0)
	{
		return false;
	}
	else j_left = int(center_x - ellipse_b2);
	if (int(center_x + ellipse_b2) > image_mat.cols - 1)
	{
		return false;
	}
	else j_right = int(center_x + ellipse_b2);

	int coed_bites_num;
	switch (code_bites_type)
	{
	case CodeBites15:
		coed_bites_num = 15;
		break;
	case CodeBites12:
		coed_bites_num = 12;
		break;
	}

	QList<QPointF> ellipse_ring_points;
	QList<int> ellipse_ring_gray_value;

	//将圆环分成15段，每段取20个点的值
	int per_num_point = 10;

	for (int i = 0; i < coed_bites_num * per_num_point; i++)
	{
		float x = cos(float(i) / coed_bites_num / per_num_point * 2 * M_PI);
		float y = sin(float(i) / coed_bites_num / per_num_point * 2 * M_PI);

		float ell_cor_r1_x = ellipse_a1 * x;
		float ell_cor_r1_y = ellipse_b1 * y;
		float ell_cor_r2_x = ellipse_a2 * x;
		float ell_cor_r2_y = ellipse_b2 * y;

		float car_cor_r1_x = ell_cor_r1_x * cos(angle_in_pi) - ell_cor_r1_y * sin(angle_in_pi) + center_x;
		float car_cor_r1_y = ell_cor_r1_x * sin(angle_in_pi) + ell_cor_r1_y * cos(angle_in_pi) + center_y;
		float car_cor_r2_x = ell_cor_r2_x * cos(angle_in_pi) - ell_cor_r2_y * sin(angle_in_pi) + center_x;
		float car_cor_r2_y = ell_cor_r2_x * sin(angle_in_pi) + ell_cor_r2_y * cos(angle_in_pi) + center_y;

		QPoint point1 = QPoint(int(car_cor_r1_x), int(car_cor_r1_y));
		QPoint point2 = QPoint(int(car_cor_r2_x), int(car_cor_r2_y));

		QList<int> gray_list = GetALineGrayList(image_mat, point1, point2);
		//int mid_value = MIdValue(gray_list);
		int mid_value = AverageOfList(gray_list);

		ellipse_ring_points.append(QPointF(x, y));
		ellipse_ring_gray_value.append(mid_value);
	}

	////找单位圆编码的边界，一维卷积模板半宽为3
	int half_length = 5;
	int edge_index = 0;
	float ans_max = 0;
	for (int i = 0; i < ellipse_ring_points.size() - 2 * half_length; i++)
	{
		float ans = 0;
		for (int j = 0; j < half_length; j++)
		{
			ans += -ellipse_ring_gray_value.at(i + j) + ellipse_ring_gray_value.at(i + j + half_length);
		}
		if (ans > ans_max)
		{
			ans_max = ans;
			edge_index = i + half_length;
		}
	}

	//阈值
	int delta_M = 20;
	if (ans_max / half_length > delta_M)
	{
		return false;
	}
	else
		return true;
}
bool ImageDetectMethod::FindSubPixelPosOfCircleCenter_opnecv(const cv::Mat& image_mat, 
	float center_x, float center_y, float ellipse_a, float ellipse_b,
	float angle_in_pi, const std::vector<cv::Point>& contour_points,
	float& sub_pixel_center_x, float& sub_pixel_center_y,
	std::vector<cv::Point2f>* subpixel_edge_points /*= NULL*/,
	MarkPointColorType color_type/* = BlackDownWhiteUp*/)
{
	float ellipse_a2 = ellipse_a * 1.5;
	float ellipse_b2 = ellipse_b * 1.5;

	int i_top = 0;
	int i_bottom = 0;
	int j_left = 0;
	int j_right = 0;
	if (int(center_y - ellipse_b2) < 0)
	{
		i_top = 0;
	}
	else i_top = int(center_y - ellipse_b2);
	if (int(center_y + ellipse_b2) > image_mat.rows - 1)
	{
		i_bottom = image_mat.rows - 1;
	}
	else i_bottom = int(center_y + ellipse_b2);
	if (int(center_x - ellipse_b2) < 0)
	{
		j_left = 0;
	}
	else j_left = int(center_x - ellipse_b2);
	if (int(center_x + ellipse_b2) > image_mat.cols - 1)
	{
		j_right = image_mat.cols - 1;
	}
	else j_right = int(center_x + ellipse_b2);

	cv::Mat crop_image = image_mat(cv::Range(i_top, i_bottom), cv::Range(j_left, j_right));

	//cv::namedWindow("canny", cv::WINDOW_FREERATIO);
	//imshow("canny", crop_image);
	////imwrite("E:\\canny.bmp",show_contours_image);
	//cv::waitKey(0);

	cv::SimpleBlobDetector::Params params;
	params.minThreshold = 15;
	params.maxThreshold = 240;
	params.thresholdStep = 5;
	cv::Ptr<cv::SimpleBlobDetector> detector = cv::SimpleBlobDetector::create(params);
	params.filterByArea = true;
	params.maxArea = 3.1415 * ellipse_a * ellipse_b / 2.0 / 2.0 * 4.0;
	params.minArea = 3.1415 * ellipse_a * ellipse_b /2.0 / 2.0 / 4.0;
	params.blobColor = 0;
	std::vector<cv::KeyPoint> keypoints;
	detector->detect(crop_image, keypoints);
	if (keypoints.size() == 0)
	{
		return false;
	}
	else
	{
		int min_index = 0;
		float min_value = 1e20;
		for (int ii = 0; ii < keypoints.size(); ii++)
		{
			if (abs((float)i_top + keypoints[ii].pt.x - center_y) + abs((float)j_left + keypoints[ii].pt.y - center_x) < min_value)
			{
				min_value = abs((float)i_top + keypoints[ii].pt.x - center_y) + abs((float)j_left + keypoints[ii].pt.y - center_x);
				min_index = ii;
			}
		}
		sub_pixel_center_y = (float)i_top + keypoints[min_index].pt.y;
		sub_pixel_center_x = (float)j_left + keypoints[min_index].pt.x;
	}
	//if (color_type == BlackDownWhiteUp || color_type == Uncertainty)
	//{

	//}


	return true;
}


QRect ImageDetectMethod::GetEllipseROIRect(const cv::Mat& image_mat, float center_x, float center_y, float ellipse_a, float ellipse_b, float angle_in_pi)
{
	int delta = 0;

	int i_top = int(center_y - ellipse_b - delta);
	if (i_top < 0)
	{
		i_top = 0;
	}

	int i_bottom = ceil(center_y + ellipse_b + delta);
	if (i_bottom > image_mat.rows - 1)
	{
		i_bottom = image_mat.rows - 1;
	}
	int j_left = int(center_x - ellipse_b - delta);
	if (j_left < 0)
	{
		j_left = 0;
	}

	int j_right = ceil(center_x + ellipse_b + delta);
	if (j_right > image_mat.cols - 1)
	{
		j_right = image_mat.cols - 1;
	}

	return QRect(j_left, i_top, j_right - j_left, i_bottom - i_top);
}




int ImageDetectMethod::AverageOfList(QList<int>& list_value)
{
	int num = 0;
	int aver = 0;
	for (int i = 0; i < list_value.size(); i++)
	{
		aver += list_value[i];
		num++;
	}
	aver /= num;
	return aver;
}

QList<int> ImageDetectMethod::GetALineGrayList(cv::Mat image_mat, QPoint point1, QPoint point2)
{
	QList<int> gray_list;

	if (point1.x() == point2.x())
	{
		int j = point1.x();
		if (point1.y() < point2.y())
		{
			for (int i = point1.y(); i < point2.y(); i++)
			{
				gray_list.append(image_mat.at<uchar>(i, j));
			}
		}
		else
		{
			for (int i = point2.y(); i < point1.y(); i++)
			{
				gray_list.append(image_mat.at<uchar>(i, j));
			}
		}
	}
	else
	{
		double k = double(point2.y() - point1.y()) / double(point2.x() - point1.x());
		if (abs(k) <= 1)
		{
			if (point1.x() < point2.x())
			{
				for (int j = point1.x(); j <= point2.x(); j++)
				{
					int i = int(k * (j - point1.x()) + point1.y());
					gray_list.append(image_mat.at<uchar>(i, j));
				}
			}
			else
			{
				// 				for (int j = point2.x();j<point1.x();j++)
				// 				{
				// 					int i = int(k*(j-point1.x())+point1.y());
				// 					gray_list.append(image_mat.at<uchar>(i,j));
				// 				}

				for (int j = point1.x(); j >= point2.x(); j--)
				{
					int i = int(k * (j - point1.x()) + point1.y());
					gray_list.append(image_mat.at<uchar>(i, j));
				}
			}
		}
		else
		{
			k = double(point2.x() - point1.x()) / double(point2.y() - point1.y());
			if (point1.y() < point2.y())
			{
				for (int i = point1.y(); i <= point2.y(); i++)
				{
					int j = int(k * (i - point1.y()) + point1.x());
					gray_list.append(image_mat.at<uchar>(i, j));
				}
			}
			else
			{
				// 				for (int i = point2.y();i<=point1.y();i++)
				// 				{
				// 					int j = int(k*(i-point1.y())+point1.x());
				// 					gray_list.append(image_mat.at<uchar>(i,j));
				// 				}
				for (int i = point1.y(); i >= point2.y(); i--)
				{
					int j = int(k * (i - point1.y()) + point1.x());
					gray_list.append(image_mat.at<uchar>(i, j));
				}
			}
		}
	}
	return gray_list;
}
double ImageDetectMethod::MeanValue(QList<float> value_list)
{
	double mean_v = 0;
	double mean_vn = 0;
	for (int ii =0;ii< value_list.size();ii++)
	{
		mean_v += value_list[ii];
		mean_vn += 1;
	}
	if (mean_vn)
	{
		mean_v /= mean_vn;
	}
	return mean_v;
}
int ImageDetectMethod::MIdValue(QList<int> value_list)
{
	int mid_value = 0;
	int n = value_list.size();
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n - 1 - i; j++)
		{
			if (value_list.value(i) > value_list.value(i + 1))
			{
				value_list.swapItemsAt(j, j + 1);
			}
		}
	}
	if (n / 2 == 0)
	{
		mid_value = (value_list.value(n / 2 - 1) + value_list.value(n / 2)) / 2;
	}
	else
		mid_value = value_list.value(n / 2);

	return mid_value;
}
int ImageDetectMethod::Change2To10(QList<int> list_code2)
{
	int ans10 = 0;
	int n = list_code2.size();
	for (int i = 0; i < n; i++)
	{
		ans10 = ans10 + int(pow(2.0, n - i - 1) * list_code2.value(i));
	}
	return ans10;
}

bool ImageDetectMethod::FindCircleGrid(cv::Mat ori_image, int h_num, int v_num,
	int h_offset, int v_offset,
	int h_mid_length, int v_mid_length,
	std::vector<cv::Point2f>& corners, std::vector<uchar>& sign_list,
	int& useful_corner_num,
	std::vector<std::vector<cv::Point2f>>& edge_points /*= std::vector<vector<cv::Point2f>>()*/,
	float max_ratio, float ratio, float min_radius/* = 5*/, float max_radius /*= 50*/, float ellipse_error_pixel /*=0.5*/, int min_arc,
	int min_points, int min_contour_num,
	DetectContoursMethod image_process_method /*= OTSU_Method*/,
	SubPixelPosMethod subpixel_pos_method /*=Gray_Centroid*/)
{
	if (!ori_image.data)        // 判断图片调入是否成功
		return false;        // 调入图片失败则退出

	int useful_corner_num_max = 0;
	std::vector<cv::Point2f> corners_max;
	cv::Mat processed_image_mat;
	std::vector<std::vector<cv::Point>> contours;
	QList<QList<float>> ellipse_pars;
	std::vector<std::vector<cv::Point>>contours_for_key;
	contours_for_key.resize(h_num * v_num);
	QList<QList<float>> ellipse_pars_ori;
	//1.图像预处理
	processed_image_mat = ori_image.clone();
	ImagePreprocess(ori_image, processed_image_mat);

	//2.边缘检测，存入闭合轮廓
	DetectClosedContours(processed_image_mat, contours, image_process_method);

	//3.轮廓筛选，尺寸，形状等准则，圆
	FilterEllipseContours(contours, min_radius, max_radius,
		ellipse_error_pixel, ellipse_pars, max_ratio, min_points, min_contour_num, min_arc);
	QList<QList<float>> ellipse_pars_copy = ellipse_pars;
	FilterEllipseContours_for_distance(contours, ellipse_pars);
	if (ellipse_pars.size() < h_num * v_num * 0.75)
	{
		return false;
	}

	FilterEllipseContoursForCSICalibrationPlane(processed_image_mat,
		ellipse_pars_copy, ellipse_pars, ellipse_pars_ori, ratio);

	std::vector<std::vector<cv::Point>> contours_copy = contours;
	contours.clear();
	for (int ii = 0; ii < contours_copy.size(); ii++)
	{
		bool need_con_c = false;
		for (int pp = 0; pp < ellipse_pars.size(); pp++)
		{
			if (ii == ellipse_pars[pp][6])
			{
				need_con_c = true;
				break;
			}
		}
		if (need_con_c)
		{
			contours.push_back(contours_copy[ii]);
		}
	}
	//cv::Mat II = cv::Mat::zeros(ori_image.size(), CV_8UC3);
	//cv::cvtColor(ori_image, II, CV_GRAY2BGR);
	////cv::drawContours(II, contours_copy, -1, cv::Scalar(0, 0, 255), 10);
	//cv::drawContours(II, contours, -1, cv::Scalar(0, 255, 0), 5);
	//cv::namedWindow("canny", cv::WINDOW_NORMAL);
	//cv::imshow("canny", II);
	//cv::waitKey(0);
	std::vector<cv::Point2f> points_temp;
	corners.clear();
	for (int ii = 0; ii < ellipse_pars.size(); ii++)
	{
		points_temp.push_back(cv::Point2f(ellipse_pars[ii][0], ellipse_pars[ii][1]));
	}

	//cv::CirclesGridFinderParameters parameters;
	//parameters.gridType = cv::CirclesGridFinderParameters::SYMMETRIC_GRID;
	//CirclesGridClusterFinder circlesGridClusterFinder(parameters);
	//circlesGridClusterFinder.findGrid(points_temp, cv::Size(h_num,v_num), corners);
	//useful_corner_num = corners.size();
	//if (useful_corner_num != h_num * v_num)
	//{
	//	bool isValid = false;
	//	const int attempts = 5000;
	//	const size_t minHomographyPoints = 4;
	//	cv::Mat H;
	//	for (int i = 0; i < attempts; i++)
	//	{
	//		corners.clear();
	//		CirclesGridFinder boxFinder(cv::Size(h_num, v_num), points_temp, parameters);
	//		try
	//		{
	//			bool isFound = boxFinder.findHoles();
	//			if (isFound)
	//			{
	//				boxFinder.getHoles(corners);
	//				isValid = true;
	//				break;  // done, return result
	//			}
	//		}
	//		catch (const cv::Exception& e)
	//		{
	//			CV_UNUSED(e);
	//			// nothing, next attempt
	//		}

	//		boxFinder.getHoles(corners);
	//		if (i != attempts - 1)
	//		{
	//			if (corners.size() < minHomographyPoints)
	//				break;
	//			H = CirclesGridFinder::rectifyGrid(boxFinder.getDetectedGridSize(), corners, points_temp, points_temp);
	//		}
	//	}

	//	if (!corners.empty() && !H.empty())  // undone rectification
	//	{
	//		cv::Mat orgPointsMat;
	//		transform(corners, orgPointsMat, H.inv());
	//		convertPointsFromHomogeneous(orgPointsMat, corners);
	//		useful_corner_num = corners.size();
	//	}
	//}
	if (ellipse_pars_ori.size() > 50)
	{
		return false;
	}
	for (int kk = 0; kk < ellipse_pars_ori.size() - 2; kk++)
	{
		for (int pp = kk + 1; pp < ellipse_pars_ori.size() - 1; pp++)
		{
			for (int qq = pp + 1; qq < ellipse_pars_ori.size(); qq++)
			{
				std::vector<int> orident_point_index_list;
				orident_point_index_list.push_back(kk);
				orident_point_index_list.push_back(pp);
				orident_point_index_list.push_back(qq);
				int original_point_index = 0, X_axis_point_index = 0, Y_axis_point_index = 0;
				float d_max = 0;
				for (int i = 0; i < orident_point_index_list.size() - 1; i++)
				{
					for (int j = i + 1; j < orident_point_index_list.size(); j++)
					{
						float x1 = ellipse_pars_ori[orident_point_index_list[i]][0];
						float y1 = ellipse_pars_ori[orident_point_index_list[i]][1];
						float x2 = ellipse_pars_ori[orident_point_index_list[j]][0];
						float y2 = ellipse_pars_ori[orident_point_index_list[j]][1];
						float length = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
						if (length > d_max)
						{
							d_max = length;
							X_axis_point_index = orident_point_index_list[i];
							Y_axis_point_index = orident_point_index_list[j];
						}
					}
				}
				for (int i = 0; i < orident_point_index_list.size(); i++)
				{
					if (orident_point_index_list[i] != X_axis_point_index && orident_point_index_list[i] != Y_axis_point_index)
					{
						original_point_index = orident_point_index_list[i];
					}
				}

				//***判断两个直角边上两方向点之间的圆点个数，确定正确的XY轴
				std::vector<float> line_X_points_length;
				std::vector<float> line_Y_points_length;
				std::vector<cv::Point2f> line_X_points;
				std::vector<cv::Point2f> line_Y_points;
				cv::Point2f original_point(ellipse_pars_ori[original_point_index][0], ellipse_pars_ori[original_point_index][1]);
				cv::Point2f X_axis_point(ellipse_pars_ori[X_axis_point_index][0], ellipse_pars_ori[X_axis_point_index][1]);
				cv::Point2f Y_axis_point(ellipse_pars_ori[Y_axis_point_index][0], ellipse_pars_ori[Y_axis_point_index][1]);
				float X_line_A, X_line_B, X_line_C, Y_line_A, Y_line_B, Y_line_C;
				LineEquation(original_point, X_axis_point, X_line_A, X_line_B, X_line_C);
				LineEquation(original_point, Y_axis_point, Y_line_A, Y_line_B, Y_line_C);

				double mean_ab = (sqrt(ellipse_pars_ori[original_point_index][2] * ellipse_pars_ori[original_point_index][3]) / 2.0
					+ sqrt(ellipse_pars_ori[X_axis_point_index][2] * ellipse_pars_ori[X_axis_point_index][3]) / 2.0
					+ sqrt(ellipse_pars_ori[Y_axis_point_index][2] * ellipse_pars_ori[Y_axis_point_index][3]) / 2.0) / 3.0;
				double dis_r_x = 0;
				double dis_r_y = 0;
				double L_x = sqrt(pow(original_point.x - X_axis_point.x, 2) + pow(original_point.y - X_axis_point.y, 2));
				double L_y = sqrt(pow(original_point.x - Y_axis_point.x, 2) + pow(original_point.y - Y_axis_point.y, 2));
				std::vector<int> Line_x_exist;
				std::vector<int> Line_y_exist;
				if (h_mid_length > v_mid_length)
				{
					if (L_x > L_y)
					{
						dis_r_x = L_x / (double)h_mid_length;
						dis_r_y = L_y / (double)v_mid_length;
					}
					else
					{
						dis_r_x = L_x / (double)v_mid_length;
						dis_r_y = L_y / (double)h_mid_length;
					}
				}
				else
				{
					if (L_x > L_y)
					{
						dis_r_x = L_x / (double)v_mid_length;
						dis_r_y = L_y / (double)h_mid_length;
					}
					else
					{
						dis_r_x = L_x / (double)h_mid_length;
						dis_r_y = L_y / (double)v_mid_length;
					}
				}

				for (int i = 0; i < ellipse_pars.size(); i++)
				{
					if (ellipse_pars[i][6]== ellipse_pars_ori[original_point_index][6])
					{
						continue;
					}
					if (ellipse_pars[i][6] == ellipse_pars_ori[X_axis_point_index][6])
					{
						continue;
					}
					if (ellipse_pars[i][6] == ellipse_pars_ori[Y_axis_point_index][6])
					{
						continue;
					}

					if (sqrt(ellipse_pars[i][2] * ellipse_pars[i][2]) / 2.0 < mean_ab * 0.5 || sqrt(ellipse_pars[i][2] * ellipse_pars[i][2]) / 2.0 > mean_ab * 2)
					{
						continue;
					}
					double L2ori = sqrt(pow(original_point.x - ellipse_pars[i][0], 2) + pow(original_point.y - ellipse_pars[i][1], 2));

					float d_max = ellipse_pars[i][2];

					cv::Point2f p(ellipse_pars[i][0], ellipse_pars[i][1]);
					float dis_to_X_line = PointTpLineDistance(p, X_line_A, X_line_B, X_line_C);
					float dis_to_Y_line = PointTpLineDistance(p, Y_line_A, Y_line_B, Y_line_C);

					float dot_in_X_line = (p.x - original_point.x) * (p.x - X_axis_point.x) + (p.y - original_point.y) * (p.y - X_axis_point.y);
					float dot_in_Y_line = (p.x - original_point.x) * (p.x - Y_axis_point.x) + (p.y - original_point.y) * (p.y - Y_axis_point.y);

					if (dis_to_X_line < d_max && dot_in_X_line < 0)
					{
						line_X_points_length.push_back(L2ori);
						line_X_points.push_back(p);
					}

					if (dis_to_Y_line < d_max && dot_in_Y_line < 0)
					{
						line_Y_points.push_back(p);
						line_Y_points_length.push_back(L2ori);
					}
				}
				////Approximately 1st re_calcu
				//std::vector<float> line_X_points_length_copy = line_X_points_length;
				//std::vector<float> line_Y_points_length_copy = line_Y_points_length;
				//std::sort(line_X_points_length.begin(), line_X_points_length.end());
				//std::sort(line_Y_points_length.begin(), line_Y_points_length.end());
				//std::vector<float> line_X_points_delta;
				//std::vector<float> line_Y_points_delta;
				//std::vector<float> line_X_points_delta2;
				//std::vector<float> line_Y_points_delta2;
				//std::vector<int> line_X_points_new_index;
				//line_X_points_new_index.resize(line_X_points_length_copy.size());
				//std::vector<int> line_Y_points_new_index;
				//line_Y_points_new_index.resize(line_Y_points_length_copy.size());
				//for (int ii = 0; ii < line_X_points_length_copy.size(); ii++)
				//{
				//	int cur_index = 0;
				//	for (int jj = 0; jj < line_X_points_length_copy.size(); jj++)
				//	{
				//		if (line_X_points_length[ii] == line_X_points_length_copy[jj])
				//		{
				//			cur_index = jj;
				//			break;
				//		}
				//	}
				//	line_X_points_new_index[ii] = cur_index;
				//}
				//for (int ii = 0; ii < line_Y_points_length_copy.size(); ii++)
				//{
				//	int cur_index = 0;
				//	for (int jj = 0; jj < line_Y_points_length_copy.size(); jj++)
				//	{
				//		if (line_Y_points_length[ii] == line_Y_points_length_copy[jj])
				//		{
				//			cur_index = jj;
				//			break;
				//		}
				//	}
				//	line_Y_points_new_index[ii] = cur_index;
				//}
				//line_X_points_delta.push_back(line_X_points_length[line_X_points_new_index[0]]);
				//for (int ii = 1; ii < line_X_points_length_copy.size(); ii++)
				//{
				//	line_X_points_delta.push_back(line_X_points_length[line_X_points_new_index[ii]] - line_X_points_length[line_X_points_new_index[ii - 1]]);
				//}
				//line_X_points_delta.push_back(sqrt(pow(X_axis_point.x - line_X_points[line_X_points_new_index[line_X_points_new_index.size() - 1]].x, 2) 
				//	+ pow(X_axis_point.y - line_X_points[line_X_points_new_index[line_X_points_new_index.size() - 1]].y, 2)));
				//line_Y_points_delta.push_back(line_Y_points_length[line_Y_points_new_index[0]]);
				//for (int ii = 1; ii < line_Y_points_length_copy.size(); ii++)
				//{
				//	line_Y_points_delta.push_back(line_Y_points_length[line_Y_points_new_index[ii]] - line_Y_points_length[line_Y_points_new_index[ii - 1]]);
				//}
				//line_Y_points_delta.push_back(sqrt(pow(Y_axis_point.x - line_Y_points[line_Y_points_new_index[line_Y_points_new_index.size() - 1]].x, 2)
				//	+ pow(Y_axis_point.y - line_Y_points[line_Y_points_new_index[line_Y_points_new_index.size() - 1]].y, 2)));

				//for (int ii = 1; ii < line_X_points_delta.size(); ii++)
				//{
				//	line_X_points_delta2.push_back(line_X_points_delta[ii] - line_X_points_delta[ii - 1]);
				//}
				//for (int ii = 1; ii < line_Y_points_delta.size(); ii++)
				//{
				//	line_Y_points_delta2.push_back(line_Y_points_delta[ii] - line_Y_points_delta[ii - 1]);
				//}
				//size_t n = line_X_points_delta2.size() / 2;
				//std::nth_element(line_X_points_delta2.begin(), line_X_points_delta2.begin() + n, line_X_points_delta2.end());
				//float X_delta = line_X_points_delta2[n];
				//n = line_Y_points_delta2.size() / 2;
				//std::nth_element(line_Y_points_delta2.begin(), line_Y_points_delta2.begin() + n, line_Y_points_delta2.end());
				//float Y_delta = line_Y_points_delta2[n];

				//n = line_X_points_delta.size() / 2;
				//std::nth_element(line_X_points_delta.begin(), line_X_points_delta.begin() + n, line_X_points_delta.end());
				//float X_Ldelta = line_X_points_delta[n];
				//n = line_Y_points_delta.size() / 2;
				//std::nth_element(line_Y_points_delta.begin(), line_Y_points_delta.begin() + n, line_Y_points_delta.end());
				//float Y_Ldelta = line_Y_points_delta[n];

				//std::vector<int> can_pass_X;
				//for (int ii = 1; ii < line_X_points_delta.size(); ii++)
				//{

				//}
				corners.clear();
				sign_list.clear();
				corners.resize(h_num * v_num);
				for (int sc = 0; sc < h_num * v_num; sc++)
				{
					corners[sc].x = 0;
					corners[sc].y = 0;
				}
				sign_list.resize(h_num * v_num);
				useful_corner_num = 0;

				if ((line_X_points.size() == h_mid_length - 1 && line_Y_points.size() == v_mid_length - 1)
					||(line_X_points.size() == v_mid_length - 1 && line_Y_points.size() == h_mid_length - 1))
				{
					if (line_X_points.size() == v_mid_length - 1 && line_Y_points.size() == h_mid_length - 1)
					{
						std::swap(X_axis_point_index, Y_axis_point_index);
						std::swap(line_X_points, line_Y_points);
						std::swap(X_axis_point, Y_axis_point);
					}

					//右手系判断，向量叉乘大于0
					cv::Point2f vector1(X_axis_point.x - original_point.x, X_axis_point.y - original_point.y);
					cv::Point2f vector2(Y_axis_point.x - original_point.x, Y_axis_point.y - original_point.y);
					if (vector1.x * vector2.y - vector1.y * vector2.x > 0)
					{
						continue;
					}
					std::vector<cv::Point2f> src_points, dst_points;
					src_points.push_back(original_point);
					src_points.push_back(X_axis_point);
					src_points.push_back(Y_axis_point);
					dst_points.push_back(cv::Point2f(h_offset, v_offset));
					dst_points.push_back(cv::Point2f(h_offset + h_mid_length, v_offset));
					dst_points.push_back(cv::Point2f(h_offset, v_offset + v_mid_length));
					for (int i = 0; i < line_X_points.size() - 1; i++)
					{
						for (int j = 0; j < line_X_points.size() - i - 1; j++)
						{
							float d1 = PointToPointDistance(line_X_points[j], original_point);
							float d2 = PointToPointDistance(line_X_points[j + 1], original_point);
							if (d1 > d2)
							{
								std::swap(line_X_points[j], line_X_points[j + 1]);
							}
						}
					}
					for (int i = 0; i < line_Y_points.size() - 1; i++)
					{
						for (int j = 0; j < line_Y_points.size() - i - 1; j++)
						{
							float d1 = PointToPointDistance(line_Y_points[j], original_point);
							float d2 = PointToPointDistance(line_Y_points[j + 1], original_point);
							if (d1 > d2)
							{
								std::swap(line_Y_points[j], line_Y_points[j + 1]);
							}
						}
					}
					for (int i = 0; i < line_X_points.size(); i++)
					{
						src_points.push_back(line_X_points[i]);
						dst_points.push_back(cv::Point2f(h_offset + i + 1, v_offset));
					}
					for (int i = 0; i < line_Y_points.size(); i++)
					{
						src_points.push_back(line_Y_points[i]);
						dst_points.push_back(cv::Point2f(h_offset, v_offset + i + 1));
					}

					//
					cv::Mat H_mat = cv::findHomography(src_points, dst_points, cv::RANSAC);
					H_mat.convertTo(H_mat, CV_32F);

					//5.4 重新按顺序排列

					std::vector<float> coners_cost;
					coners_cost.resize(h_num * v_num);
					for (int ii = 0; ii < coners_cost.size(); ii++)
					{
						coners_cost[ii] = -1;
					}
					for (int n = 0; n < ellipse_pars.size(); n++)
					{
						//if (sqrt(ellipse_pars[n][2] * ellipse_pars[n][2]) / 2.0 < mean_ab * 0.5 || sqrt(ellipse_pars[n][2] * ellipse_pars[n][2]) / 2.0 > mean_ab * 2)
						//{
						//	continue;
						//}
						cv::Mat X = cv::Mat(3, 1, CV_32F);
						X.at<float>(0, 0) = ellipse_pars[n][0];
						X.at<float>(1, 0) = ellipse_pars[n][1];
						X.at<float>(2, 0) = 1;
						cv::Mat A = H_mat * X;
						float new_x = A.at<float>(0, 0) / A.at<float>(2, 0);
						float new_y = A.at<float>(1, 0) / A.at<float>(2, 0);

						int i = floor(new_y + 0.5);
						int j = floor(new_x + 0.5);
						float delta_x = abs(new_x - j);
						float delta_y = abs(new_y - i);

						if (i<1 || i>v_num || j<1 || j>h_num)
						{
							continue;
						}
						if ((i - 1)==9 && (j - 1) == 12)
						{
							int scacas = 1;
						}
						if (coners_cost[(i - 1) * h_num + j - 1] == -1)
						{
							coners_cost[(i - 1) * h_num + j - 1] = delta_x + delta_y;
							float sub_pixel_x, sub_pixel_y;
							std::vector<cv::Point2f> edge_contour;
							if (FindSubPixelPosOfCircleCenter_opnecv(processed_image_mat, ellipse_pars[n][0], ellipse_pars[n][1], ellipse_pars[n][2],
								ellipse_pars[n][3], ellipse_pars[n][4], contours_copy[ellipse_pars[n][6]], sub_pixel_x, sub_pixel_y,
								&edge_contour, Uncertainty))
							{
								ellipse_pars[n][0] = sub_pixel_x;
								ellipse_pars[n][1] = sub_pixel_y;
								corners[(i - 1) * h_num + j - 1] = cv::Point2f(sub_pixel_x, sub_pixel_y);
								sign_list[(i - 1) * h_num + j - 1] = 2;
							}
							else
							{
								corners[(i - 1) * h_num + j - 1] = cv::Point2f(ellipse_pars[n][0], ellipse_pars[n][1]);
								sign_list[(i - 1) * h_num + j - 1] = 1;
							}
							contours_for_key[(i - 1) * h_num + j - 1] = contours_copy[ellipse_pars[n][6]];
							sign_list[(i - 1) * h_num + j - 1] = 2;
							useful_corner_num++;
						}
						else if (coners_cost[(i - 1) * h_num + j - 1] < (delta_x + delta_y))
						{
							continue;
						}
						else
						{
							coners_cost[(i - 1) * h_num + j - 1] = delta_x + delta_y;
							float sub_pixel_x, sub_pixel_y;
							std::vector<cv::Point2f> edge_contour;
							if (FindSubPixelPosOfCircleCenter_opnecv(processed_image_mat, ellipse_pars[n][0], ellipse_pars[n][1], ellipse_pars[n][2],
								ellipse_pars[n][3], ellipse_pars[n][4], contours_copy[ellipse_pars[n][6]], sub_pixel_x, sub_pixel_y,
								&edge_contour, Uncertainty))
							{
								ellipse_pars[n][0] = sub_pixel_x;
								ellipse_pars[n][1] = sub_pixel_y;
								corners[(i - 1) * h_num + j - 1] = cv::Point2f(sub_pixel_x, sub_pixel_y);
								sign_list[(i - 1) * h_num + j - 1] = 2;
							}
							else
							{
								corners[(i - 1) * h_num + j - 1] = cv::Point2f(ellipse_pars[n][0], ellipse_pars[n][1]);
								sign_list[(i - 1) * h_num + j - 1] = 1;
							}
							contours_for_key[(i - 1) * h_num + j - 1] = contours_copy[ellipse_pars[n][6]];
						}
					}
					if (useful_corner_num> useful_corner_num_max)
					{
						useful_corner_num_max = useful_corner_num;
						corners_max = corners;
					}
					if (useful_corner_num == (h_num * v_num))
					{
						//cv::Mat II = cv::Mat::zeros(ori_image.size(), CV_8UC3);
						//cv::cvtColor(ori_image, II, CV_GRAY2BGR);
						////cv::drawContours(II, contours_copy, -1, cv::Scalar(0, 0, 255), 10);
						//cv::drawContours(II, contours_for_key, -1, cv::Scalar(0, 255, 0), 5);
						//cv::namedWindow("cannys", cv::WINDOW_NORMAL);
						//cv::imshow("cannys", II);
						//cv::waitKey(0);
						return true;
					}
				}
			}
		}
	}

	if (useful_corner_num == (h_num * v_num))
	{
		return true;
	}
	else if (useful_corner_num_max >= (h_num * v_num) * 0.75)
	{
		useful_corner_num = useful_corner_num_max;
		corners = corners_max;
		return true;
	}
	return false;
}

bool ImageDetectMethod::FilterEllipseContours_for_distance(std::vector<std::vector<cv::Point>> contours, QList<QList<float>>& ellipse_pars)
{
	std::vector<int> remove_index;
	for (int i = 0; i < ellipse_pars.size() - 1; i++)
	{
		bool need_con = false;
		for (int pp = 0; pp < remove_index.size(); pp++)
		{
			if (i == remove_index[pp])
			{
				need_con = true;
				break;
			}
		}
		if (need_con)
		{
			continue;
		}
		for (int j = i + 1; j < ellipse_pars.size(); j++)
		{
			bool need_con_c = false;
			for (int pp = 0; pp < remove_index.size(); pp++)
			{
				if (j == remove_index[pp])
				{
					need_con_c = true;
					break;
				}
			}
			if (need_con_c)
			{
				continue;
			}
			float x1 = ellipse_pars[i][0];
			float y1 = ellipse_pars[i][1];
			float x2 = ellipse_pars[j][0];
			float y2 = ellipse_pars[j][1];
			float length_of_2_points = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));

			if (length_of_2_points < std::min(ellipse_pars[i][3], ellipse_pars[j][3]))
			{
				if (contours[ellipse_pars[i][6]].size() > contours[ellipse_pars[j][6]].size())
				{
					remove_index.push_back(j);
				}
				else
				{
					remove_index.push_back(i);
				}

				//if (ellipse_pars[i][3] > ellipse_pars[j][3])
				//{
				//	remove_index.push_back(j);
				//}
				//else
				//{
				//	remove_index.push_back(i);
				//}
			}
		}
	}
	QList<QList<float>> ellipse_pars_copy = ellipse_pars;
	ellipse_pars.clear();
	for (int ii = 0; ii < ellipse_pars_copy.size(); ii++)
	{
		bool need_con_c = false;
		for (int pp = 0; pp < remove_index.size(); pp++)
		{
			if (ii == remove_index[pp])
			{
				need_con_c = true;
				break;
			}
		}
		if (!need_con_c)
		{
			ellipse_pars.push_back(ellipse_pars_copy[ii]);
		}
	}
	return true;
}
bool ImageDetectMethod::FilterEllipseContoursForCSICalibrationPlane(const cv::Mat& image_mat,
	QList<QList<float>>& ellipse_pars_all,
	QList<QList<float>>& ellipse_pars, QList<QList<float>>& ellipse_pars_ori, float ratio_k)
{
	QList<float> gray_value_std_list;
	for (int i = 0; i < ellipse_pars.size(); i++)
	{
		if (EllipseGrayJudgeForPointCSI_is2Circle(ellipse_pars_all, ellipse_pars[i], ratio_k))
		{
			ellipse_pars_ori.append(ellipse_pars[i]);
		}
	}
	if (ellipse_pars_ori.size() < 3)
	{
		return false;
	}
	return true;
}

void ImageDetectMethod::LineEquation(cv::Point2f p1, cv::Point2f p2, float& A, float& B, float& C)
{
	float x1 = p1.x;
	float y1 = p1.y;
	float x2 = p2.x;
	float y2 = p2.y;

	if (x2 == x1)
	{
		A = 1;
		B = 0;
		C = -x1;

		return;
	}
	else
	{
		float k = (y2 - y1) / (x2 - x1);
		float b = y1 - k * x1;

		A = k;
		B = -1;
		C = b;
		return;
	}
}

float ImageDetectMethod::PointTpLineDistance(cv::Point2f p, float A, float B, float C)
{
	float d = abs(A * p.x + B * p.y + C) / sqrt(A * A + B * B);
	return d;
}

float ImageDetectMethod::PointToPointDistance(cv::Point2f p1, cv::Point2f p2)
{
	return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
}

void writeCSV(std::string filename, cv::Mat m)
{
	std::ofstream myfile;
	myfile.open(filename.c_str());
	myfile << cv::format(m, cv::Formatter::FMT_CSV) << std::endl;
	myfile.close();
}
//bool ImageDetectMethod::detectCheckerboard(cv::Mat I, std::vector<cv::Point2d>& corner_points, double sigma, double peakThreshold, bool highDistortion, bool usePartial)
//{
//	if (I.channels() == 3)
//	{
//		cv::cvtColor(I, I, CV_BGR2GRAY);
//	}
//	if (I.type() != CV_32F)
//	{
//		if (I.type() == CV_8U)
//		{
//			I.convertTo(I, CV_32F, 1 / 255.0);
//		}
//		else
//		{
//			I.convertTo(I, CV_32F);
//		}
//	}
//	cv::Mat cxy, c45, Ix, Iy, Ixy, I_45_45;
//	secondDerivCornerMetric(I, cxy, c45, Ix, Iy, Ixy, I_45_45, sigma, highDistortion);
//
//	cv::Mat Ix2 = Ix.mul(Ix);
//	cv::Mat Iy2 = Iy.mul(Iy);
//	cv::Mat Ixy2 = Ix.mul(Iy);
//	cv::GaussianBlur(Ix2, Ix2, cv::Size(7, 7), 1.5);
//	cv::GaussianBlur(Iy2, Iy2, cv::Size(7, 7), 1.5);
//	cv::GaussianBlur(Ixy, Ix2, cv::Size(7, 7), 1.5);
//	return true;
//}

bool ImageDetectMethod::detectCheckerboard(cv::Mat I, std::vector<cv::Point2f>& corner_points, cv::Size& wh_check, double peakThreshold, bool highDistortion, bool usePartial)
{
	if (I.channels() == 3)
	{
		cv::cvtColor(I, I, CV_BGR2GRAY);
	}
	coder::array<double, 2U> imagePoints;
	coder::array<unsigned char, 2U> b_I;
	b_I.set_size(I.rows, I.cols);
	for (int idx0{ 0 }; idx0 < b_I.size(0); idx0++) {
		for (int idx1{ 0 }; idx1 < b_I.size(1); idx1++) {
			b_I[idx0 + b_I.size(0) * idx1] = I.at<uchar>(idx0, idx1);
		}
	}
	double boardSize[2];
	boolean_T B1_tmp;
	boolean_T imagesUsed;
	get_chessborad_pixel(b_I, peakThreshold, highDistortion, usePartial, imagePoints,
		boardSize, &imagesUsed);
	if (imagesUsed)
	{
		corner_points.clear();
		for (int idx0{ 0 }; idx0 < imagePoints.size(0); idx0++)
		{
			corner_points.push_back(cv::Point2d(imagePoints[idx0], imagePoints[idx0 + imagePoints.size(0)]));
			wh_check.width = boardSize[0];
			wh_check.height = boardSize[1];
		}
		return true;
	}
	else
	{
		return false;
	}
}

bool ImageDetectMethod::ImageDetectMethod::detectcodecircle_graysearch(cv::Mat ori_image, cv::Mat& code_point_mat,
	float ratio_k, float ratio_k1, float ratio_k2,
	float min_radius, float max_radius, float ellipse_error_pixel /*=0.5*/, 
	MarkPointColorType color_type /*=BlackDownWhiteUp*/, CodePointBitesType code_bites_type /*=CodeBites15*/
	, SubPixelPosMethod subpixel_pos_method
	, int search_gray_min, int search_gray_max, int search_gray_step
	, double max_aspect_ratio, int min_points)
{
	if (!ori_image.data)
		return false;
	if (ori_image.channels() == 3)
	{
		cv::cvtColor(ori_image, ori_image, CV_BGR2GRAY);
	}
	cv::Mat processed_image_mat;
	QList<QList<float>> ellipse_pars_all;

	//图像预处理
	ImagePreprocess(ori_image, processed_image_mat);//高斯滤波前处理，C
	if (color_type == MarkPointColorType::BlackDownWhiteUp || color_type == MarkPointColorType::Uncertainty)
	{
		for (int gray_threshold = search_gray_min; gray_threshold < search_gray_max; gray_threshold = gray_threshold + search_gray_step)
		{
			QList<QList<float>> ellipse_pars;
			std::vector<cv::Vec4i> hierarchy;
			std::vector<std::vector<cv::Point>> contours;
			cv::Mat threshold_image;
			cv::threshold(processed_image_mat, threshold_image, gray_threshold, 255, cv::THRESH_BINARY);
			findContours(threshold_image, contours, hierarchy, cv::RETR_LIST, cv::CHAIN_APPROX_NONE, cv::Point(0, 0));

			//轮廓筛选，尺寸，形状等准则
			FilterEllipseContours(contours, min_radius, max_radius,
				ellipse_error_pixel, ellipse_pars, max_aspect_ratio, min_points);

			//进一步筛选 ，用于编码点和非编码
			FilterEllipseContoursForCodePoint(processed_image_mat, ratio_k, ratio_k1, ratio_k2,
				ellipse_pars);

			int* default_id_array_ptr = NULL;
			int default_id_array_size;
			default_id_array_ptr = ReturnDefualtIdArray(default_id_array_size, code_bites_type);
			if (default_id_array_ptr == NULL)
			{
				return false;
			}

			int uncodePoint_id = 0;
			for (int i = 0; i < ellipse_pars.size(); i++)
			{
				int code_id;
				bool is_decode = Decoding20140210(processed_image_mat, code_id, ellipse_pars[i][0], ellipse_pars[i][1], ellipse_pars[i][2], ellipse_pars[i][3], ellipse_pars[i][4],
					ratio_k1, ratio_k2, BlackDownWhiteUp, code_bites_type);

				bool is_code_point = false;
				if (is_decode == true)
				{
					for (int j = 0; j < default_id_array_size; j++)
					{
						int id = *(default_id_array_ptr + j);
						if (code_id == *(default_id_array_ptr + j))
						{
							is_code_point = true;

							ellipse_pars[i].append(j);
							ellipse_pars[i].append(1);

							break;
						}
					}
				}

				if (is_code_point == false)
				{
					bool 	is_uncodepoint = true;
					if (is_uncodepoint == true)
					{
						ellipse_pars[i].append(uncodePoint_id);
						ellipse_pars[i].append(0);
						uncodePoint_id++;
					}
					else
					{
						ellipse_pars.removeAt(i);
						i--;
					}
				}
			}
			std::vector<std::vector<cv::Point2f>> subpixel_edge_contours;
			for (int i = 0; i < ellipse_pars.size(); i++)
			{
				float sub_pixel_x, sub_pixel_y;
				std::vector<cv::Point2f> edge_contour;
				if (FindSubPixelPosOfCircleCenter_opnecv(processed_image_mat, ellipse_pars[i][0], ellipse_pars[i][1], ellipse_pars[i][2],
					ellipse_pars[i][3], ellipse_pars[i][4], contours[ellipse_pars[i][6]], sub_pixel_x, sub_pixel_y,
					&edge_contour, Uncertainty))
				{
					ellipse_pars[i][0] = sub_pixel_x;
					ellipse_pars[i][1] = sub_pixel_y;
				}
				subpixel_edge_contours.push_back(edge_contour);
			}
			for (int i = 0; i < ellipse_pars.size(); i++)
			{
				ellipse_pars_all.append(ellipse_pars[i]);
			}
		}
	}
	if (color_type == MarkPointColorType::WhiteDownBlackUp || color_type == MarkPointColorType::Uncertainty)
	{
		for (int gray_threshold = search_gray_min; gray_threshold < search_gray_max; gray_threshold = gray_threshold + search_gray_step)
		{
			QList<QList<float>> ellipse_pars;
			std::vector<cv::Vec4i> hierarchy;
			std::vector<std::vector<cv::Point>> contours;
			cv::Mat threshold_image;
			cv::threshold(processed_image_mat, threshold_image, gray_threshold, 255, cv::THRESH_BINARY_INV);
			findContours(threshold_image, contours, hierarchy, cv::RETR_LIST, cv::CHAIN_APPROX_NONE, cv::Point(0, 0));

			//轮廓筛选，尺寸，形状等准则
			FilterEllipseContours(contours, min_radius, max_radius,
				ellipse_error_pixel, ellipse_pars, max_aspect_ratio, min_points);

			//进一步筛选 ，用于编码点和非编码
			FilterEllipseContoursForCodePoint(processed_image_mat, ratio_k, ratio_k1, ratio_k2,
				ellipse_pars);

			int* default_id_array_ptr = NULL;
			int default_id_array_size;
			default_id_array_ptr = ReturnDefualtIdArray(default_id_array_size, code_bites_type);
			if (default_id_array_ptr == NULL)
			{
				return false;
			}

			int uncodePoint_id = 0;
			for (int i = 0; i < ellipse_pars.size(); i++)
			{
				int code_id;
				bool is_decode = Decoding20140210(processed_image_mat, code_id, ellipse_pars[i][0], ellipse_pars[i][1], ellipse_pars[i][2], ellipse_pars[i][3], ellipse_pars[i][4],
					ratio_k1, ratio_k2, WhiteDownBlackUp, code_bites_type);

				bool is_code_point = false;
				if (is_decode == true)
				{
					for (int j = 0; j < default_id_array_size; j++)
					{
						int id = *(default_id_array_ptr + j);
						if (code_id == *(default_id_array_ptr + j))
						{
							is_code_point = true;

							ellipse_pars[i].append(j);
							ellipse_pars[i].append(1);

							break;
						}
					}
				}

				if (is_code_point == false)
				{
					bool 	is_uncodepoint = true;
					if (is_uncodepoint == true)
					{
						ellipse_pars[i].append(uncodePoint_id);
						ellipse_pars[i].append(0);
						uncodePoint_id++;
					}
					else
					{
						ellipse_pars.removeAt(i);
						i--;
					}
				}
			}
			std::vector<std::vector<cv::Point2f>> subpixel_edge_contours;
			for (int i = 0; i < ellipse_pars.size(); i++)
			{
				float sub_pixel_x, sub_pixel_y;
				std::vector<cv::Point2f> edge_contour;
				if (FindSubPixelPosOfCircleCenter_opnecv(processed_image_mat, ellipse_pars[i][0], ellipse_pars[i][1], ellipse_pars[i][2],
					ellipse_pars[i][3], ellipse_pars[i][4], contours[ellipse_pars[i][6]], sub_pixel_x, sub_pixel_y,
					&edge_contour, Uncertainty))
				{
					ellipse_pars[i][0] = sub_pixel_x;
					ellipse_pars[i][1] = sub_pixel_y;
				}
				subpixel_edge_contours.push_back(edge_contour);
			}
			for (int i = 0; i < ellipse_pars.size(); i++)
			{
				ellipse_pars_all.append(ellipse_pars[i]);
			}
		}
	}

	code_point_mat = cv::Mat();
	for (int i = 0; i < ellipse_pars_all.size(); i++)
	{
		float a[7] = { ellipse_pars_all[i][7],ellipse_pars_all[i][0],ellipse_pars_all[i][1],ellipse_pars_all[i][5],
			ellipse_pars_all[i][2],ellipse_pars_all[i][3],ellipse_pars_all[i][4] };
		cv::Mat mat = cv::Mat(1, 7, CV_32F, a);
		if (ellipse_pars_all[i][8] > 0)
		{
			code_point_mat.push_back(mat);
		}
	}
	return true;
}
bool ImageDetectMethod::detectcodecircle(cv::Mat ori_image, cv::Mat& code_point_mat, std::vector<std::vector<cv::Point>>& contours_pixel,
	float ratio_k, float ratio_k1, float ratio_k2,
	float min_radius, float max_radius, float ellipse_error_pixel /*=0.5*/,
	MarkPointColorType color_type /*= BlackDownWhiteUp*/, CodePointBitesType code_bites_type /*=CodeBites15*/,
	DetectContoursMethod image_process_method /*= CANNY_Method*/,
	SubPixelPosMethod subpixel_pos_method,
	double max_aspect_ratio, int min_points, int min_contour_num,
	float delta_Mt, float fore_stdDev, float back_stdDev)
{
	if (!ori_image.data)        // 判断图片调入是否成功
		return false;        // 调入图片失败则退出

	if (ori_image.channels() == 3)
	{
		cv::cvtColor(ori_image, ori_image, CV_BGR2GRAY);
	}
	contours_pixel.clear();
	QList<QList<float>> ellipse_pars_all;
	cv::Mat processed_image_mat;
	ImagePreprocess(ori_image, processed_image_mat);

	if (color_type == MarkPointColorType::BlackDownWhiteUp || color_type == MarkPointColorType::Uncertainty)
	{
		std::vector<std::vector<cv::Point>> contours;
		QList<QList<float>> ellipse_pars;
		DetectClosedContours(processed_image_mat, contours, image_process_method);
		FilterEllipseContours(contours, min_radius, max_radius,
			ellipse_error_pixel, ellipse_pars, max_aspect_ratio, min_points, min_contour_num);
		FilterEllipseContoursForCodePoint(processed_image_mat, ratio_k, ratio_k1, ratio_k2,
			ellipse_pars, delta_Mt, fore_stdDev, back_stdDev);
		int* default_id_array_ptr = NULL;
		int default_id_array_size;
		default_id_array_ptr = ReturnDefualtIdArray(default_id_array_size, code_bites_type);
		if (default_id_array_ptr == NULL)
		{
			return false;
		}

		int uncodePoint_id = 0;
		for (int i = 0; i < ellipse_pars.size(); i++)
		{
			int code_id;
			bool is_decode = Decoding20140210(processed_image_mat, code_id, ellipse_pars[i][0], ellipse_pars[i][1], ellipse_pars[i][2], ellipse_pars[i][3], ellipse_pars[i][4],
				ratio_k1, ratio_k2, BlackDownWhiteUp, code_bites_type, delta_Mt);

			bool is_code_point = false;
			if (is_decode == true)
			{
				for (int j = 0; j < default_id_array_size; j++)
				{
					int id = *(default_id_array_ptr + j);
					if (code_id == *(default_id_array_ptr + j))
					{
						is_code_point = true;

						ellipse_pars[i].append(j);
						ellipse_pars[i].append(1);

						break;
					}
				}
			}

			if (is_code_point == false)
			{
				bool 	is_uncodepoint = true;
				if (is_uncodepoint == true)
				{
					ellipse_pars[i].append(uncodePoint_id);
					ellipse_pars[i].append(0);
					uncodePoint_id++;
				}
				else
				{
					ellipse_pars.removeAt(i);
					i--;
				}
			}
		}

		std::vector<std::vector<cv::Point2f>> subpixel_edge_contours;
		for (int i = 0; i < ellipse_pars.size(); i++)
		{
			float sub_pixel_x, sub_pixel_y;
			std::vector<cv::Point2f> edge_contour;
			if (FindSubPixelPosOfCircleCenter_opnecv(processed_image_mat, ellipse_pars[i][0], ellipse_pars[i][1], ellipse_pars[i][2],
				ellipse_pars[i][3], ellipse_pars[i][4], contours[ellipse_pars[i][6]], sub_pixel_x, sub_pixel_y,
				&edge_contour, Uncertainty))
			{
				ellipse_pars[i][0] = sub_pixel_x;
				ellipse_pars[i][1] = sub_pixel_y;
			}
			subpixel_edge_contours.push_back(edge_contour);
		}

		for (int i = 0; i < ellipse_pars.size(); i++)
		{
			if (ellipse_pars[i][8] > 0)
			{
				ellipse_pars_all.append(ellipse_pars[i]);
				contours_pixel.push_back(contours[ellipse_pars[i][6]]);
			}
		}
	}
	if (color_type == MarkPointColorType::WhiteDownBlackUp || color_type == MarkPointColorType::Uncertainty)
	{
		std::vector<std::vector<cv::Point>> contours;
		QList<QList<float>> ellipse_pars;
		processed_image_mat = 255 - processed_image_mat;
		DetectClosedContours(processed_image_mat, contours, image_process_method);
		FilterEllipseContours(contours, min_radius, max_radius,
			ellipse_error_pixel, ellipse_pars, max_aspect_ratio, min_points, min_contour_num);
		FilterEllipseContoursForCodePoint(processed_image_mat, ratio_k, ratio_k1, ratio_k2,
			ellipse_pars, delta_Mt, fore_stdDev, back_stdDev);
		int* default_id_array_ptr = NULL;
		int default_id_array_size;
		default_id_array_ptr = ReturnDefualtIdArray(default_id_array_size, code_bites_type);
		if (default_id_array_ptr == NULL)
		{
			return false;
		}

		int uncodePoint_id = 0;
		for (int i = 0; i < ellipse_pars.size(); i++)
		{
			int code_id;
			bool is_decode = Decoding20140210(processed_image_mat, code_id, ellipse_pars[i][0], ellipse_pars[i][1], ellipse_pars[i][2], ellipse_pars[i][3], ellipse_pars[i][4],
				ratio_k1, ratio_k2, WhiteDownBlackUp, code_bites_type, delta_Mt);

			bool is_code_point = false;
			if (is_decode == true)
			{
				for (int j = 0; j < default_id_array_size; j++)
				{
					int id = *(default_id_array_ptr + j);
					if (code_id == *(default_id_array_ptr + j))
					{
						is_code_point = true;

						ellipse_pars[i].append(j);
						ellipse_pars[i].append(1);

						break;
					}
				}
			}

			if (is_code_point == false)
			{
				bool 	is_uncodepoint = true;
				if (is_uncodepoint == true)
				{
					ellipse_pars[i].append(uncodePoint_id);
					ellipse_pars[i].append(0);
					uncodePoint_id++;
				}
				else
				{
					ellipse_pars.removeAt(i);
					i--;
				}
			}
		}

		std::vector<std::vector<cv::Point2f>> subpixel_edge_contours;
		for (int i = 0; i < ellipse_pars.size(); i++)
		{
			float sub_pixel_x, sub_pixel_y;
			std::vector<cv::Point2f> edge_contour;
			if (FindSubPixelPosOfCircleCenter_opnecv(processed_image_mat, ellipse_pars[i][0], ellipse_pars[i][1], ellipse_pars[i][2],
				ellipse_pars[i][3], ellipse_pars[i][4], contours[ellipse_pars[i][6]], sub_pixel_x, sub_pixel_y,
				&edge_contour, Uncertainty))
			{
				ellipse_pars[i][0] = sub_pixel_x;
				ellipse_pars[i][1] = sub_pixel_y;
			}
			subpixel_edge_contours.push_back(edge_contour);
		}

		for (int i = 0; i < ellipse_pars.size(); i++)
		{
			if (ellipse_pars[i][8] > 0)
			{
				ellipse_pars_all.append(ellipse_pars[i]);
				contours_pixel.push_back(contours[ellipse_pars[i][6]]);
			}
		}
	}
//ellipse_pars - n*6  center_x,center_y,r_a,r_b,angle_inPI,ellipse_error,contours_index,ID,code_type(0- uncode,1- code)
	code_point_mat = cv::Mat();
	for (int i = 0; i < ellipse_pars_all.size(); i++)
	{
		float a[7] = { ellipse_pars_all[i][7],ellipse_pars_all[i][0],ellipse_pars_all[i][1],ellipse_pars_all[i][5],
			ellipse_pars_all[i][2],ellipse_pars_all[i][3],ellipse_pars_all[i][4] };
		cv::Mat mat = cv::Mat(1, 7, CV_32F, a);
		code_point_mat.push_back(mat);
	}
	return true;
}