#ifndef IMAGEDETECTMETHOD_H
#define IMAGEDETECTMETHOD_H

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

	//*********************编码点和非编码点检测********************************************
	//image_file_name - 输入图像路径名
	//ratio_k - 输入灰度准则的比例，非编码点外径与内径之比
	//ratio_k1 - 输入编码环内径与中间定位圆的内径之比
	//ratio_k2 - 输入编码环外径与中间定位圆的内径之比
	//min_radius - 圆最小半径
	//max_radius - 圆最大半径
	//color_type - 标志点白底黑点还是黑底白点
	//code_bites_type - 编码点编码位数， 一般默认常用 15_8 15段编码8个1
	//code_point_mat - 输出编码点信息矩阵  n*7 , id,x,y,quality, r_a,r_b,angle_in_pi
	//uncode_point_mat - 输出编码点信息矩阵 n*7 , id,x,y,quality, r_a,r_b,angle_in_pi
	static bool CodeAndUncodePointDetect(QString image_file_name, cv::Mat& code_point_mat, cv::Mat& uncode_point_mat,
		float ratio_k, float ratio_k1, float ratio_k2,
		float min_radius, float max_radius, float ellipse_error_pixel = 0.5,
		MarkPointColorType color_type = BlackDownWhiteUp, CodePointBitesType code_bites_type = CodeBites15,
		DetectContoursMethod image_process_method = OTSU_Method, SubPixelPosMethod subpixel_pos_method = Gray_Centroid);

	//*******************圆心定位精度测试用**********************
	static bool TestCircleCenterPosAccuracy(QString image_file_name, cv::Mat& code_point_mat, cv::Mat& uncode_point_mat,
		float ratio_k, float ratio_k1, float ratio_k2,
		float min_radius, float max_radius, float ellipse_error_pixel = 0.5,
		MarkPointColorType color_type = BlackDownWhiteUp, CodePointBitesType code_bites_type = CodeBites15,
		DetectContoursMethod image_process_method = OTSU_Method, SubPixelPosMethod subpixel_pos_method = Gray_Centroid);

	//标定板角点检测
	//corners-输出角点坐标，sign_list-输出点是否有效
	//棋盘
	static bool FindChessGridOpencv(QString image_file_name, int h_num, int v_num, std::vector<cv::Point2f>& corners, std::vector<uchar>& sign_list, int& useful_corner_num);

	//CSI标定板
	static bool FindCircleGrid(cv::Mat ori_image, int h_num, int v_num, int h_offset, int v_offset, int h_mid_length,
		int v_mid_length, std::vector<cv::Point2f>& corners, std::vector<uchar>& sign_list, int& useful_corner_num,
		std::vector<std::vector<cv::Point2f>>& edge_points,
		float max_ratio = 3, float ratio = 0.5,
		float min_radius = 5, float max_radius = 50, float ellipse_error_pixel = 0.5 , int min_arc = 240, int min_points = 6, int min_contour_num = 8,
		DetectContoursMethod image_process_method = CANNY_Method, SubPixelPosMethod subpixel_pos_method = Gauss_Curve_Fit);


	//编码点自标定，特征点检测
	static bool FindCodePointsForSelfCalibrationByHartly(QStringList image_file_name_list, std::vector<std::vector<cv::Point2f>>& code_points_list,
		float ratio_k = 2, float ratio_k1 = 2.4, float ratio_k2 = 4,
		float min_radius = 5, float max_radius = 50, float ellipse_error_pixel = 0.5,
		MarkPointColorType color_type = BlackDownWhiteUp, CodePointBitesType code_bites_type = CodeBites15,
		DetectContoursMethod image_process_method = OTSU_Method, SubPixelPosMethod subpixel_pos_method = Gray_Centroid);
	
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
	static cv::Mat Adaptive_Canny(const cv::Mat& src, int apertureSize, bool L2gradient);
	static cv::Mat _get_canny_result(const cv::Mat& map);
	static void _hysteresis_thresholding(std::deque<int>& mapIndicesX, std::deque<int>& mapIndicesY, cv::Mat& map);
	static void _non_maximum_suppression(const cv::Mat& NMSImage, cv::Mat& map, std::deque<int>& mapIndicesX,
		std::deque<int>& mapIndicesY, int low, int high);
	static void _calculate_hysteresis_threshold_value(const cv::Mat& dx, const cv::Mat& dy, const cv::Mat& magnitudes,
		const cv::Mat& angles, cv::Mat& NMSImage, int& low, int& high);
	static void _sobel_gradient(const cv::Mat& mat, cv::Mat& dx, cv::Mat& dy, cv::Mat& magnitudes, cv::Mat& angles,
		int apertureSize, bool L2gradient);
	//检测闭合的轮廓线并存入序列
	static bool DetectClosedContours(const cv::Mat& image_mat, std::vector<std::vector<cv::Point>>& contours, DetectContoursMethod image_process_method);
	//二值图，闭合边缘追踪
	static bool TraceEdge(int start_y, int start_x, int y, int x,
		std::vector<cv::Point>* edge_vector, cv::Mat* cannyed_image_mat, cv::Mat* result_mat, bool* is_contour_closed);
	//自己写的Canny算子
	static bool SelfCannyMethod(const cv::Mat& ori_image_mat, cv::Mat& output_image_mat, float canny_low_thresh, float canny_high_thresh,
		std::vector<std::vector<cv::Point>>& edge_point_list);
	//canny算子边缘追踪，高低阈值
	static bool FindCannyEdge(int start_y, int start_x, int y, int x,
		cv::Mat* low_threshold_mat, cv::Mat* high_threshold_mat, cv::Mat* search_sign_mat, std::vector<cv::Point>* edge_point_vector, bool* is_contour_closed);

	//条件准则约束，筛选边缘轮廓
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

	static bool EllipseGrayJudgeForPointCSI2(const cv::Mat& image_mat, float center_x, float center_y,
		float ellipse_a, float ellipse_b, float angle_in_pi, float ratio_k);
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

	static bool FindSubPixelPosOfCircleCenter20140210(const cv::Mat& image_mat, float center_x, float center_y, float ellipse_a, float ellipse_b,
		float angle_in_pi, const std::vector<cv::Point>& contour_points,
		float& sub_pixel_center_x, float& sub_pixel_center_y, std::vector<cv::Point2f>* subpixel_edge_points = NULL, SubPixelPosMethod subPixel_method = NoSubPixel_Match, MarkPointColorType color_type = BlackDownWhiteUp);

	//获取椭圆局部感兴趣区域
	static QRect GetEllipseROIRect(const cv::Mat& image_mat, float center_x, float center_y, float ellipse_a, float ellipse_b, float angle_in_pi);

	//**********Otsu方法图像二值化，返回二值化阈值,支持高比特
	static float DICOTSU20140215(const cv::Mat& ori_image, int total_level = 256);

	//******如果没有指定，黑底白点或白底黑点时，判断是哪种类型
	static MarkPointColorType JudgeTargetColorType(const cv::Mat& sub_mat, float center_x_insubmat, float center_y_insubmat, float ellipse_a, float ellipse_b, float angle_in_pi);

	//******一种获得局部阈值的放法， 阈值 = (目标区域平均灰度 + 背景区域平局灰度)/2
	static float CalThresholdInSubsetMat(const cv::Mat& sub_mat, float center_x_insubmat, float center_y_insubmat, float ellipse_a, float ellipse_b, float angle_in_pi);

	//*****一种获取局部阈值的方法，轮廓边缘灰度平均值
	static float CalThresholdInSubsetMat2(const cv::Mat& image_mat, const std::vector<cv::Point>& contour_points);

	//*********灰度重心法计算局部重心
	static bool CalCentriodBySubsetMat(const cv::Mat& subset_mat, float thresh_old, float& sub_center_x, float& sub_center_y,
		MarkPointColorType color_type = BlackDownWhiteUp, SubPixelPosMethod subPixel_method = Gray_Centroid);

	//**********获得重心法灰度权重
	static float CalWeightOfCentriod(float gray_value, float gray_threshold, MarkPointColorType color_type = BlackDownWhiteUp,
		SubPixelPosMethod subPixel_method = Gray_Centroid);

	//*******计算点的灰度梯度幅值和角度（G）
	static void CalGrayGradientBySobel(const cv::Mat& image_mat, int i, int j, float* gray_gradient, float* gradient_sita = NULL);

	//*********曲率滤波，计算连续边缘点的曲率
	static void CalCurvatureFromEdgePoints(const std::vector<cv::Point2f>& edge_points, std::vector<float>& curvature_vector);

	//************根据拟合误差，筛选部分点重新选点进行拟合
	static void ReduceBadEllipseFitPoints(std::vector<cv::Point2f>& edge_points, float center_x, float center_y, float ellipse_a, float ellipse_b, float angle_in_pi);

	//************测试用显示轮廓
	static void ShowContous(int image_width, int image_height, std::vector<std::vector<cv::Point>>& contours, QList<QList<float>>& ellipse_pars);

	//**************平均值****
	static int AverageOfList(QList<int>& list_value);

	//*****************canny算子自适应阈值Otsu法确定高阈值Th2，Th1= 0.5Th2；
	//返回处于哪一等级时类间方差最大
	static double OTSUForCannyAdaptiveHighTh(cv::Mat ori_image, int total_level = 256);

	//*********************获取由椭圆圆心发出射线与编码环两个交点之间的线段上的灰度值******************
	//image_mat-输入原始图像
	//point1-输入交点1坐标
	//point2-输入交点2坐标
	//返回一系列数
	static QList<int> GetALineGrayList(cv::Mat image_mat, QPoint point1, QPoint point2);

	static double MeanValue(QList<float> value_list);
	//*************************获取一系列数的中值*********************************
	//value_list - 输入数组
	//返回中值
	static int MIdValue(QList<int> value_list);

	//***********二进制转十进制，输入list_code2二进制数，输出return十进制值**********************
	static int Change2To10(QList<int> list_code2);

	//
	static void WirteMatToFile(cv::Mat mat, QString file_name = "mat");

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

	//static bool secondDerivCornerMetric(cv::Mat I, cv::Mat& cxy, cv::Mat& c45, cv::Mat& Ix, cv::Mat& Iy, cv::Mat& Ixy, cv::Mat& I_45_45, double sigma = 2.0, bool highDistortionl = false);
	
	//static bool find_corner_peaks(cv::Mat metric, cv::Mat& loc, double peakThreshold = 0.15);
};

#endif // IMAGEDETECTMETHOD_H
