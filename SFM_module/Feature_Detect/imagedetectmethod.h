#ifndef IMAGEDETECTMETHOD_H
#define IMAGEDETECTMETHOD_H


/********************
Originally written by Professor Dong Shuai, with extensive detailed modifications by Yin Zhuoyi.
*******************/
#include <QObject>
#include "CMMHelp.h"
#include "CodeID.h"
#include<opencv2/opencv.hpp>
#include<vector>
#include <fstream>
#include "core/core.hpp"
#include <algorithm>
#include<ceres/ceres.h>
#include <ceres/problem.h>
#include "Matlab2C/rtwtypes.h"
#include <cstddef>
#include <cstdlib>
#include "Matlab2C/imedge_2d.h"
#include "Matlab2C/get_chessborad_pixel.h"
#include "Matlab2C/get_chessborad_pixel_terminate.h"
#include "Matlab2C/rt_nonfinite.h"
#include "Matlab2C/coder_array.h"
#include "Opencv/circlesgrid.hpp"

//#include "edge/imedge_2d.h"
//#include "edge/imedge_2d_terminate.h"

class Fit_ellipse_ceres
{
public:
	Fit_ellipse_ceres(const Eigen::Vector2d img_pixel) :
		m_img_pixel(img_pixel) {}

	template <typename T>
	bool operator()(const T* coff, T* residuals) const
	{
		T tr_x = (m_img_pixel.x() - coff[0]) * cos(coff[4]) + (m_img_pixel.y() - coff[1]) * sin(coff[4]);
		T tr_y = -(m_img_pixel.x() - coff[0]) * sin(coff[4]) + (m_img_pixel.y() - coff[1]) * cos(coff[4]);

		T alfa = atan2(tr_y, tr_x);

		T r = coff[2] * coff[3] / sqrt(coff[2] * coff[2] * sin(alfa) * sin(alfa) + coff[3] * coff[3] * cos(alfa) * cos(alfa));
		residuals[0] = sqrt(tr_x * tr_x + tr_y * tr_y) - r;

		return true;
	}

	// 生成误差函数
	static ceres::CostFunction* Create(Eigen::Vector2d img_pixel)
	{
		return (new ceres::AutoDiffCostFunction<Fit_ellipse_ceres, 1, 5>
			(new Fit_ellipse_ceres(img_pixel)));
	}
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
	const Eigen::Vector2d m_img_pixel;
};

class ImageDetectMethod : public QObject
{
	Q_OBJECT

public:
	ImageDetectMethod(QObject* parent);
	~ImageDetectMethod();

	//CSI标定板
	static bool FindCircleGrid(cv::Mat ori_image, int h_num, int v_num, int h_offset, int v_offset, int h_mid_length,
		int v_mid_length, std::vector<cv::Point2f>& corners, std::vector<uchar>& sign_list, int& useful_corner_num,
		std::vector<std::vector<cv::Point2f>>& edge_points,
		float max_ratio = 3, float ratio = 0.5,
		float min_radius = 5, float max_radius = 50, float ellipse_error_pixel = 0.5 , int min_arc = 240, int min_points = 6, int min_contour_num = 8,
		DetectContoursMethod image_process_method = CANNY_Method, SubPixelPosMethod subpixel_pos_method = Gauss_Curve_Fit);

	//棋盘标定板角点识别
	static bool detectCheckerboard(cv::Mat I, std::vector<cv::Point2f>& corner_points, cv::Size& wh_check, double peakThreshold = 0.15, bool highDistortion = false, bool usePartial = false);
	static bool detectcodecircle(cv::Mat ori_image, cv::Mat& code_point_mat, std::vector<std::vector<cv::Point>> & contours_pixel,
		float ratio_k, float ratio_k1, float ratio_k2,
		float min_radius, float max_radius, float ellipse_error_pixel = 0.5,
		MarkPointColorType color_type = BlackDownWhiteUp, CodePointBitesType code_bites_type = CodeBites15,
		DetectContoursMethod image_process_method = OTSU_Method, SubPixelPosMethod subpixel_pos_method = Gray_Centroid
		, double max_aspect_ratio = 2, int min_points = 6, int min_contour_num = 8,
		float delta_Mt = 50, float fore_stdDev = 100, float back_stdDev = 100);
	static bool detectcodecircle_graysearch(cv::Mat ori_image, cv::Mat& code_point_mat,
		float ratio_k, float ratio_k1, float ratio_k2,
		float min_radius, float max_radius, float ellipse_error_pixel /*=0.5*/, MarkPointColorType color_type = BlackDownWhiteUp, 
		CodePointBitesType code_bites_type =CodeBites15, SubPixelPosMethod subpixel_pos_method = Gray_Centroid
		, int search_gray_min = 20, int search_gray_max = 235, int search_gray_step = 20
		, double max_aspect_ratio = 2, int min_points = 6);

private:
	//图像预处理
	static bool ImagePreprocess(const cv::Mat& ori_image_mat, cv::Mat& processed_image_mat);
	static bool DetectClosedContours(const cv::Mat& image_mat, std::vector<std::vector<cv::Point>>& contours, DetectContoursMethod image_process_method);
	
	//ellipse_pars - n*6  center_x,center_y,r_a,r_b,angle_inPI,ellipse_error,contours_index,ID
	static bool FilterEllipseContours(const std::vector<std::vector<cv::Point>>& contours, int min_radius_pixel, int max_radius_pixel, float ellipse_error_pixel,
		QList<QList<float>>& ellipse_pars, double max_aspect_ratio = 2, int min_points = 6, int min_contour_num = 8
		, int min_arc = 240);

	static float Contours_Length(std::vector<cv::Point> C);
	static float Contours_arc(std::vector<cv::Point> C, cv::Point2f center);

	static float ErrorDROfEllipseFit(float center_x, float center_y, float ellipse_a, float ellipse_b, float ellipse_angle,
		float x, float  y);
	static float LeastSquareErrorOfEllipseFit(float center_x, float center_y, float ellipse_a, float ellipse_b, float ellipse_angle,
		float x, float  y);

	//进一步筛选出适用于编码点解码的圆点
	static bool FilterEllipseContoursForCodePoint(const cv::Mat& image_mat, float ratio_k, float ratio_k1, float ratio_k2,
		QList<QList<float>>& ellipse_pars,
		float delta_Mt = 50, float fore_stdDev = 100, float back_stdDev = 100);

	static bool EllipseGrayJudgeForPointCSI(const cv::Mat& image_mat, float center_x, float center_y,
		float ellipse_a, float ellipse_b, float angle_in_pi, float ratio_k,
		float& out_std);
	static bool EllipseGrayJudgeForPointCSI_is2Circle(
		QList<QList<float>> ellipse_pars_all, QList<float> ellipse_pars_now, float rati_k);

	static bool EllipseGrayJudgeForCodePoint(const cv::Mat& image_mat, float center_x, float center_y,
		float ellipse_a, float ellipse_b, float angle_in_pi, float ratio_k,
		float& out_foreground_mean, float& out_background_mean, float& out_foreground_stdDev, float& out_background_stdDev,
		float delta_Mt = 50, float fore_stdDev = 100, float back_stdDev = 100);

	//获取编码值数组的指针，
	static int* ReturnDefualtIdArray(int& array_size, CodePointBitesType code_bites_type = CodeBites15);

	//解码,编码点解码
	static bool Decoding20140210(const cv::Mat& image_mat, int& out_put_code_id, float center_x, float center_y, float ellipse_a, float ellipse_b, float angle_in_pi,
		float ratio_k1 = 2.4, float ratio_k2 = 4, MarkPointColorType color_type = BlackDownWhiteUp, CodePointBitesType code_bites_type = CodeBites15,
		double thread_value = 0);

	static bool CalculateRealCodeID20140210(QList<int> in_put_code_list, QList<int>& out_put_code_list, int& out_put_code_ID);
	//非编码点判断
	static bool UncodePointCheck(const cv::Mat& image_mat, float center_x, float center_y, float ellipse_a, float ellipse_b, float angle_in_pi,
		float ratio_k = 2, MarkPointColorType color_type = BlackDownWhiteUp, CodePointBitesType code_bites_type = CodeBites15);

	//亚像素定位
	static bool FindSubPixelPosOfCircleCenter_opnecv(const cv::Mat& image_mat, float center_x, float center_y, float ellipse_a, float ellipse_b,
		float angle_in_pi, const std::vector<cv::Point>& contour_points,
		float& sub_pixel_center_x, float& sub_pixel_center_y, std::vector<cv::Point2f>* subpixel_edge_points = NULL, MarkPointColorType color_type = BlackDownWhiteUp);

	//获取椭圆局部感兴趣区域
	static QRect GetEllipseROIRect(const cv::Mat& image_mat, float center_x, float center_y, float ellipse_a, float ellipse_b, float angle_in_pi);


	//******如果没有指定，黑底白点或白底黑点时，判断是哪种类型
	static MarkPointColorType JudgeTargetColorType(const cv::Mat& sub_mat, float center_x_insubmat, float center_y_insubmat, float ellipse_a, float ellipse_b, float angle_in_pi);



	//*********************获取由椭圆圆心发出射线与编码环两个交点之间的线段上的灰度值******************
	//image_mat-输入原始图像
	//point1-输入交点1坐标
	//point2-输入交点2坐标
	//返回一系列数
	static QList<int> GetALineGrayList(cv::Mat image_mat, QPoint point1, QPoint point2);

	//**************平均值****
	static int AverageOfList(QList<int>& list_value);

	static double MeanValue(QList<float> value_list);
	//*************************获取一系列数的中值*********************************
	//value_list - 输入数组
	//返回中值
	static int MIdValue(QList<int> value_list);

	//***********二进制转十进制，输入list_code2二进制数，输出return十进制值**********************
	static int Change2To10(QList<int> list_code2);

	static bool FilterEllipseContours_for_distance(std::vector<std::vector<cv::Point>> contours, QList<QList<float>>& ellipse_pars);

	//筛选出CSI类圆标定板的 圆点方向点（外圆点）+ 普通点， vector<int> orident_point_index,3个方向点的index
	static bool FilterEllipseContoursForCSICalibrationPlane(const cv::Mat& image_mat,
		QList<QList<float>>& ellipse_pars_all,
		QList<QList<float>>& ellipse_pars, QList<QList<float>>& ellipse_pars_ori, float ratio_k = 0.5);

	//计算直线方程Ax+By+C=0
	static void LineEquation(cv::Point2f p1, cv::Point2f p2, float& A, float& B, float& C);

	//计算点p到直线Ax+By+C=0距离d
	static float PointTpLineDistance(cv::Point2f p, float A, float B, float C);

	//计算点到点之间距离
	static float PointToPointDistance(cv::Point2f p1, cv::Point2f p2);

};

#endif // IMAGEDETECTMETHOD_H
