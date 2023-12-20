#pragma once
#include<ceres/ceres.h>
#include <ceres/problem.h>
#include "Eigen/Eigen"
#include<opencv2/opencv.hpp>
#include <vector>
#include "Format_Document.h"


class ProjectErrorCostFunction_PNP
{
public:
	ProjectErrorCostFunction_PNP(Eigen::Vector3d points, Eigen::Vector2d observed, Eigen::Matrix<double, 3, 3> camera_params
		, Eigen::Matrix<double, 1, 6> camera_params_dis_k, Eigen::Matrix<double, 1, 2> camera_params_dis_p, Eigen::Matrix<double, 1, 4> camera_params_dis_t) :
		m_points(points), m_observed(observed), m_camera_params(camera_params), m_camera_params_dis_k(camera_params_dis_k)
		, m_camera_params_dis_p(camera_params_dis_p), m_camera_params_dis_t(camera_params_dis_t)
	{}
	template <typename T>
	bool operator()(const T* const camQvec, const T* const camTvec, T* residuals) const
	{
		Eigen::Matrix<T, 3, 1> obj;
		obj << T(m_points(0)), T(m_points(1)), T(m_points(2));

		Eigen::Quaternion<T> qvec(camQvec);
		Eigen::Matrix<T, 3, 1> tvec;
		tvec << camTvec[0], camTvec[1], camTvec[2];

		Eigen::Matrix<T, 3, 1> obj_cam_coor = qvec.toRotationMatrix() * obj + tvec;

		T a = obj_cam_coor(0) / obj_cam_coor(2);
		T b = obj_cam_coor(1) / obj_cam_coor(2);

		T r2 = (a * a + b * b);
		T r4 = r2 * r2;
		T r6 = r2 * r4;

		T xd = a * (1.0 + T(m_camera_params_dis_k(0)) * r2 + T(m_camera_params_dis_k(1)) * r4 + T(m_camera_params_dis_k(2)) * r6)
			/ (1.0 + T(m_camera_params_dis_k(3)) * r2 + T(m_camera_params_dis_k(4)) * r4 + T(m_camera_params_dis_k(5)) * r6) +
			2.0 * T(m_camera_params_dis_p(0)) * a * b + T(m_camera_params_dis_p(1)) * (r2 + 2.0 * a * a) + T(m_camera_params_dis_t(0)) * r2 + T(m_camera_params_dis_t(2)) * r4;

		T yd = b * (1.0 + T(m_camera_params_dis_k(0)) * r2 + T(m_camera_params_dis_k(1)) * r4 + T(m_camera_params_dis_k(2)) * r6)
			/ (1.0 + T(m_camera_params_dis_k(3)) * r2 + T(m_camera_params_dis_k(4)) * r4 + T(m_camera_params_dis_k(5)) * r6) +
			2.0 * T(m_camera_params_dis_p(1)) * a * b + T(m_camera_params_dis_p(0)) * (r2 + 2.0 * b * b) + T(m_camera_params_dis_t(1)) * r2 + T(m_camera_params_dis_t(3)) * r4;

		T ud = T(m_camera_params(0,0)) * xd + T(m_camera_params(0, 1)) * yd + T(m_camera_params(0, 2));
		T vd = T(m_camera_params(1, 1)) * yd + T(m_camera_params(1, 2));

		// 残差 二维
		residuals[0] = ud - T(m_observed(0));
		residuals[1] = vd - T(m_observed(1));
		return true;
	}

	// 生成误差函数
	static ceres::CostFunction* Create(Eigen::Vector3d m_points, Eigen::Vector2d m_observed, Eigen::Matrix<double, 3, 3> m_camera_params
		, Eigen::Matrix<double, 1, 6> m_camera_params_dis_k, Eigen::Matrix<double, 1, 2> m_camera_params_dis_p, Eigen::Matrix<double, 1, 4> m_camera_params_dis_t)
	{
		return (new ceres::AutoDiffCostFunction<ProjectErrorCostFunction_PNP, 2, 4, 3>
			(new ProjectErrorCostFunction_PNP(m_points, m_observed, m_camera_params, m_camera_params_dis_k, m_camera_params_dis_p, m_camera_params_dis_t)));
	}
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:

	const Eigen::Vector3d m_points;
	const Eigen::Vector2d m_observed;
	const Eigen::Matrix<double, 3, 3> m_camera_params;
	const Eigen::Matrix<double, 1, 6> m_camera_params_dis_k;
	const Eigen::Matrix<double, 1, 2> m_camera_params_dis_p;
	const Eigen::Matrix<double, 1, 4> m_camera_params_dis_t;

};

class ProjectErrorCostFunctionCircle_ori
{
public:
	ProjectErrorCostFunctionCircle_ori(const Eigen::Vector2d Cirlce_contour, const Eigen::Matrix<double,3,4> m_camera_params) :
		m_cirlce_contour(Cirlce_contour), m_camera_params(m_camera_params) {}

	/*circle_N法向方向，自由度2
	circle_params 圆心空间坐标(3) 半径(1) 
	*/
	template <typename T>
	bool operator()(const T* circle_params, T* residuals) const
	{
		T N1 = sin(circle_params[4]);
		T N2 = cos(circle_params[4]) * sin(circle_params[5]);
		T N3 = cos(circle_params[4]) * cos(circle_params[5]);
		T x_w= (m_camera_params(0, 1) * m_camera_params(1, 3) * N3 - m_camera_params(0, 2) * m_camera_params(1, 3) * N2 - m_camera_params(0, 3) * m_camera_params(1, 1) * N3 + m_camera_params(0, 3) * m_camera_params(1, 2) * N2 + m_camera_params(1, 1) * m_camera_params(2, 3) * N3 * m_cirlce_contour.x() - m_camera_params(1, 2) * m_camera_params(2, 3) * N2 * m_cirlce_contour.x() - m_camera_params(1, 3) * m_camera_params(2, 1) * N3 * m_cirlce_contour.x() + m_camera_params(1, 3) * m_camera_params(2, 2) * N2 * m_cirlce_contour.x() + m_camera_params(0, 1) * m_camera_params(1, 2) * N1 * circle_params[0] - m_camera_params(0, 2) * m_camera_params(1, 1) * N1 * circle_params[0] - m_camera_params(0, 1) * m_camera_params(2, 3) * N3 * m_cirlce_contour.y() + m_camera_params(0, 2) * m_camera_params(2, 3) * N2 * m_cirlce_contour.y() + m_camera_params(0, 3) * m_camera_params(2, 1) * N3 * m_cirlce_contour.y() - m_camera_params(0, 3) * m_camera_params(2, 2) * N2 * m_cirlce_contour.y() + m_camera_params(0, 1) * m_camera_params(1, 2) * N2 * circle_params[1] - m_camera_params(0, 2) * m_camera_params(1, 1) * N2 * circle_params[1] + m_camera_params(0, 1) * m_camera_params(1, 2) * N3 * circle_params[2] - m_camera_params(0, 2) * m_camera_params(1, 1) * N3 * circle_params[2] + m_camera_params(1, 1) * m_camera_params(2, 2) * N1 * m_cirlce_contour.x() * circle_params[0] - m_camera_params(1, 2) * m_camera_params(2, 1) * N1 * m_cirlce_contour.x() * circle_params[0] - m_camera_params(0, 1) * m_camera_params(2, 2) * N1 * circle_params[0] * m_cirlce_contour.y() + m_camera_params(0, 2) * m_camera_params(2, 1) * N1 * circle_params[0] * m_cirlce_contour.y() + m_camera_params(1, 1) * m_camera_params(2, 2) * N2 * m_cirlce_contour.x() * circle_params[1] - m_camera_params(1, 2) * m_camera_params(2, 1) * N2 * m_cirlce_contour.x() * circle_params[1] - m_camera_params(0, 1) * m_camera_params(2, 2) * N2 * m_cirlce_contour.y() * circle_params[1] + m_camera_params(0, 2) * m_camera_params(2, 1) * N2 * m_cirlce_contour.y() * circle_params[1] + m_camera_params(1, 1) * m_camera_params(2, 2) * N3 * m_cirlce_contour.x() * circle_params[2] - m_camera_params(1, 2) * m_camera_params(2, 1) * N3 * m_cirlce_contour.x() * circle_params[2] - m_camera_params(0, 1) * m_camera_params(2, 2) * N3 * m_cirlce_contour.y() * circle_params[2] + m_camera_params(0, 2) * m_camera_params(2, 1) * N3 * m_cirlce_contour.y() * circle_params[2]) / (m_camera_params(0, 0) * m_camera_params(1, 1) * N3 - m_camera_params(0, 0) * m_camera_params(1, 2) * N2 - m_camera_params(0, 1) * m_camera_params(1, 0) * N3 + m_camera_params(0, 1) * m_camera_params(1, 2) * N1 + m_camera_params(0, 2) * m_camera_params(1, 0) * N2 - m_camera_params(0, 2) * m_camera_params(1, 1) * N1 + m_camera_params(1, 0) * m_camera_params(2, 1) * N3 * m_cirlce_contour.x() - m_camera_params(1, 0) * m_camera_params(2, 2) * N2 * m_cirlce_contour.x() - m_camera_params(1, 1) * m_camera_params(2, 0) * N3 * m_cirlce_contour.x() + m_camera_params(1, 1) * m_camera_params(2, 2) * N1 * m_cirlce_contour.x() + m_camera_params(1, 2) * m_camera_params(2, 0) * N2 * m_cirlce_contour.x() - m_camera_params(1, 2) * m_camera_params(2, 1) * N1 * m_cirlce_contour.x() - m_camera_params(0, 0) * m_camera_params(2, 1) * N3 * m_cirlce_contour.y() + m_camera_params(0, 0) * m_camera_params(2, 2) * N2 * m_cirlce_contour.y() + m_camera_params(0, 1) * m_camera_params(2, 0) * N3 * m_cirlce_contour.y() - m_camera_params(0, 1) * m_camera_params(2, 2) * N1 * m_cirlce_contour.y() - m_camera_params(0, 2) * m_camera_params(2, 0) * N2 * m_cirlce_contour.y() + m_camera_params(0, 2) * m_camera_params(2, 1) * N1 * m_cirlce_contour.y());
		T y_w = -(m_camera_params(0, 0) * m_camera_params(1, 3) * N3 - m_camera_params(0, 2) * m_camera_params(1, 3) * N1 - m_camera_params(0, 3) * m_camera_params(1, 0) * N3 + m_camera_params(0, 3) * m_camera_params(1, 2) * N1 + m_camera_params(1, 0) * m_camera_params(2, 3) * N3 * m_cirlce_contour.x() - m_camera_params(1, 2) * m_camera_params(2, 3) * N1 * m_cirlce_contour.x() - m_camera_params(1, 3) * m_camera_params(2, 0) * N3 * m_cirlce_contour.x() + m_camera_params(1, 3) * m_camera_params(2, 2) * N1 * m_cirlce_contour.x() + m_camera_params(0, 0) * m_camera_params(1, 2) * N1 * circle_params[0] - m_camera_params(0, 2) * m_camera_params(1, 0) * N1 * circle_params[0] - m_camera_params(0, 0) * m_camera_params(2, 3) * N3 * m_cirlce_contour.y() + m_camera_params(0, 2) * m_camera_params(2, 3) * N1 * m_cirlce_contour.y() + m_camera_params(0, 3) * m_camera_params(2, 0) * N3 * m_cirlce_contour.y() - m_camera_params(0, 3) * m_camera_params(2, 2) * N1 * m_cirlce_contour.y() + m_camera_params(0, 0) * m_camera_params(1, 2) * N2 * circle_params[1] - m_camera_params(0, 2) * m_camera_params(1, 0) * N2 * circle_params[1] + m_camera_params(0, 0) * m_camera_params(1, 2) * N3 * circle_params[2] - m_camera_params(0, 2) * m_camera_params(1, 0) * N3 * circle_params[2] + m_camera_params(1, 0) * m_camera_params(2, 2) * N1 * m_cirlce_contour.x() * circle_params[0] - m_camera_params(1, 2) * m_camera_params(2, 0) * N1 * m_cirlce_contour.x() * circle_params[0] - m_camera_params(0, 0) * m_camera_params(2, 2) * N1 * circle_params[0] * m_cirlce_contour.y() + m_camera_params(0, 2) * m_camera_params(2, 0) * N1 * circle_params[0] * m_cirlce_contour.y() + m_camera_params(1, 0) * m_camera_params(2, 2) * N2 * m_cirlce_contour.x() * circle_params[1] - m_camera_params(1, 2) * m_camera_params(2, 0) * N2 * m_cirlce_contour.x() * circle_params[1] - m_camera_params(0, 0) * m_camera_params(2, 2) * N2 * m_cirlce_contour.y() * circle_params[1] + m_camera_params(0, 2) * m_camera_params(2, 0) * N2 * m_cirlce_contour.y() * circle_params[1] + m_camera_params(1, 0) * m_camera_params(2, 2) * N3 * m_cirlce_contour.x() * circle_params[2] - m_camera_params(1, 2) * m_camera_params(2, 0) * N3 * m_cirlce_contour.x() * circle_params[2] - m_camera_params(0, 0) * m_camera_params(2, 2) * N3 * m_cirlce_contour.y() * circle_params[2] + m_camera_params(0, 2) * m_camera_params(2, 0) * N3 * m_cirlce_contour.y() * circle_params[2]) / (m_camera_params(0, 0) * m_camera_params(1, 1) * N3 - m_camera_params(0, 0) * m_camera_params(1, 2) * N2 - m_camera_params(0, 1) * m_camera_params(1, 0) * N3 + m_camera_params(0, 1) * m_camera_params(1, 2) * N1 + m_camera_params(0, 2) * m_camera_params(1, 0) * N2 - m_camera_params(0, 2) * m_camera_params(1, 1) * N1 + m_camera_params(1, 0) * m_camera_params(2, 1) * N3 * m_cirlce_contour.x() - m_camera_params(1, 0) * m_camera_params(2, 2) * N2 * m_cirlce_contour.x() - m_camera_params(1, 1) * m_camera_params(2, 0) * N3 * m_cirlce_contour.x() + m_camera_params(1, 1) * m_camera_params(2, 2) * N1 * m_cirlce_contour.x() + m_camera_params(1, 2) * m_camera_params(2, 0) * N2 * m_cirlce_contour.x() - m_camera_params(1, 2) * m_camera_params(2, 1) * N1 * m_cirlce_contour.x() - m_camera_params(0, 0) * m_camera_params(2, 1) * N3 * m_cirlce_contour.y() + m_camera_params(0, 0) * m_camera_params(2, 2) * N2 * m_cirlce_contour.y() + m_camera_params(0, 1) * m_camera_params(2, 0) * N3 * m_cirlce_contour.y() - m_camera_params(0, 1) * m_camera_params(2, 2) * N1 * m_cirlce_contour.y() - m_camera_params(0, 2) * m_camera_params(2, 0) * N2 * m_cirlce_contour.y() + m_camera_params(0, 2) * m_camera_params(2, 1) * N1 * m_cirlce_contour.y());
		T z_w = (m_camera_params(0, 0) * m_camera_params(1, 3) * N2 - m_camera_params(0, 1) * m_camera_params(1, 3) * N1 - m_camera_params(0, 3) * m_camera_params(1, 0) * N2 + m_camera_params(0, 3) * m_camera_params(1, 1) * N1 + m_camera_params(1, 0) * m_camera_params(2, 3) * N2 * m_cirlce_contour.x() - m_camera_params(1, 1) * m_camera_params(2, 3) * N1 * m_cirlce_contour.x() - m_camera_params(1, 3) * m_camera_params(2, 0) * N2 * m_cirlce_contour.x() + m_camera_params(1, 3) * m_camera_params(2, 1) * N1 * m_cirlce_contour.x() + m_camera_params(0, 0) * m_camera_params(1, 1) * N1 * circle_params[0] - m_camera_params(0, 1) * m_camera_params(1, 0) * N1 * circle_params[0] - m_camera_params(0, 0) * m_camera_params(2, 3) * N2 * m_cirlce_contour.y() + m_camera_params(0, 1) * m_camera_params(2, 3) * N1 * m_cirlce_contour.y() + m_camera_params(0, 3) * m_camera_params(2, 0) * N2 * m_cirlce_contour.y() - m_camera_params(0, 3) * m_camera_params(2, 1) * N1 * m_cirlce_contour.y() + m_camera_params(0, 0) * m_camera_params(1, 1) * N2 * circle_params[1] - m_camera_params(0, 1) * m_camera_params(1, 0) * N2 * circle_params[1] + m_camera_params(0, 0) * m_camera_params(1, 1) * N3 * circle_params[2] - m_camera_params(0, 1) * m_camera_params(1, 0) * N3 * circle_params[2] + m_camera_params(1, 0) * m_camera_params(2, 1) * N1 * m_cirlce_contour.x() * circle_params[0] - m_camera_params(1, 1) * m_camera_params(2, 0) * N1 * m_cirlce_contour.x() * circle_params[0] - m_camera_params(0, 0) * m_camera_params(2, 1) * N1 * circle_params[0] * m_cirlce_contour.y() + m_camera_params(0, 1) * m_camera_params(2, 0) * N1 * circle_params[0] * m_cirlce_contour.y() + m_camera_params(1, 0) * m_camera_params(2, 1) * N2 * m_cirlce_contour.x() * circle_params[1] - m_camera_params(1, 1) * m_camera_params(2, 0) * N2 * m_cirlce_contour.x() * circle_params[1] - m_camera_params(0, 0) * m_camera_params(2, 1) * N2 * m_cirlce_contour.y() * circle_params[1] + m_camera_params(0, 1) * m_camera_params(2, 0) * N2 * m_cirlce_contour.y() * circle_params[1] + m_camera_params(1, 0) * m_camera_params(2, 1) * N3 * m_cirlce_contour.x() * circle_params[2] - m_camera_params(1, 1) * m_camera_params(2, 0) * N3 * m_cirlce_contour.x() * circle_params[2] - m_camera_params(0, 0) * m_camera_params(2, 1) * N3 * m_cirlce_contour.y() * circle_params[2] + m_camera_params(0, 1) * m_camera_params(2, 0) * N3 * m_cirlce_contour.y() * circle_params[2]) / (m_camera_params(0, 0) * m_camera_params(1, 1) * N3 - m_camera_params(0, 0) * m_camera_params(1, 2) * N2 - m_camera_params(0, 1) * m_camera_params(1, 0) * N3 + m_camera_params(0, 1) * m_camera_params(1, 2) * N1 + m_camera_params(0, 2) * m_camera_params(1, 0) * N2 - m_camera_params(0, 2) * m_camera_params(1, 1) * N1 + m_camera_params(1, 0) * m_camera_params(2, 1) * N3 * m_cirlce_contour.x() - m_camera_params(1, 0) * m_camera_params(2, 2) * N2 * m_cirlce_contour.x() - m_camera_params(1, 1) * m_camera_params(2, 0) * N3 * m_cirlce_contour.x() + m_camera_params(1, 1) * m_camera_params(2, 2) * N1 * m_cirlce_contour.x() + m_camera_params(1, 2) * m_camera_params(2, 0) * N2 * m_cirlce_contour.x() - m_camera_params(1, 2) * m_camera_params(2, 1) * N1 * m_cirlce_contour.x() - m_camera_params(0, 0) * m_camera_params(2, 1) * N3 * m_cirlce_contour.y() + m_camera_params(0, 0) * m_camera_params(2, 2) * N2 * m_cirlce_contour.y() + m_camera_params(0, 1) * m_camera_params(2, 0) * N3 * m_cirlce_contour.y() - m_camera_params(0, 1) * m_camera_params(2, 2) * N1 * m_cirlce_contour.y() - m_camera_params(0, 2) * m_camera_params(2, 0) * N2 * m_cirlce_contour.y() + m_camera_params(0, 2) * m_camera_params(2, 1) * N1 * m_cirlce_contour.y());
		residuals[0] = sqrt(pow(x_w- circle_params[0],2)+ pow(y_w - circle_params[1], 2)+ pow(z_w - circle_params[2], 2))- circle_params[3];
		return true;
	}

	// 生成误差函数
	static ceres::CostFunction* Create(Eigen::Vector2d m_cirlce_contour, Eigen::Matrix<double, 3, 4> m_camera_contour)
	{
		return (new ceres::AutoDiffCostFunction<ProjectErrorCostFunctionCircle_ori, 1, 6>
			(new ProjectErrorCostFunctionCircle_ori(m_cirlce_contour, m_camera_contour)));
	}
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
	
	const Eigen::Vector2d m_cirlce_contour;
	const Eigen::Matrix<double, 3, 4> m_camera_params;
	
};
class ProjectErrorCostFunctionPinehole_SFM
{
public:
	ProjectErrorCostFunctionPinehole_SFM(const Eigen::Vector2d img_pixel) :
		m_img_pixel(img_pixel) {}

	template <typename T>
	bool operator()(const T* obj_3d, const T* const camK_F, const T* const camK_C, const T* const camK_shear,
		const T* const camDK_1, const T* const camDK_2, const T* const camDK_3, 
		const T* const camDK_4, const T* const camDK_5, const T* const camDK_6,
		const T* const camDP, const T* const camDT_2, const T* const camDT_4,
		const T* const camQvec, const T* const camTvec, T* residuals) const
	{
		// 上面的参数类型为数组

		Eigen::Matrix<T, 3, 1> obj;
		obj << obj_3d[0], obj_3d[1], obj_3d[2];

		Eigen::Quaternion<T> qvec(camQvec);
		Eigen::Matrix<T, 3, 1> tvec;
		tvec << camTvec[0], camTvec[1], camTvec[2];

		Eigen::Matrix<T, 2, 1> img = m_img_pixel.cast<T>();

		Eigen::Matrix<T, 3, 1> obj_cam_coor = qvec.toRotationMatrix() * obj + tvec;

		T a = obj_cam_coor(0) / obj_cam_coor(2);
		T b = obj_cam_coor(1) / obj_cam_coor(2);

		T r2 = (a * a + b * b);
		T r4 = r2 * r2;
		T r6 = r2 * r4;

		T xd = a * (1.0 + camDK_1[0] * r2 + camDK_2[0] * r4 + camDK_3[0] * r6)
			/ (1.0 + camDK_4[0] * r2 + camDK_5[0] * r4 + camDK_6[0] * r6) +
			2.0 * camDP[0] * a * b + camDP[1] * (r2 + 2.0 * a * a) + camDT_2[0] * r2 + +camDT_4[0] * r4;
		T yd = b * (1.0 + camDK_1[0] * r2 + camDK_2[0] * r4 + camDK_3[0] * r6)
			/ (1.0 + camDK_4[0] * r2 + camDK_5[0] * r4 + camDK_6[0] * r6) +
			2.0 * camDP[1] * a * b + camDP[0] * (r2 + 2.0 * b * b) + camDT_2[1] * r2 + +camDT_4[1] * r4;

		T ud = camK_F[0] * xd + camK_shear[0] * yd + camK_C[0];
		T vd = camK_F[1] * yd + camK_C[1];

		// 残差 二维
		residuals[0] = ud - img(0);
		residuals[1] = vd - img(1);
		
		return true;
	}

	// 生成误差函数
	static ceres::CostFunction* Create(Eigen::Vector2d img_pixel)
	{
		return (new ceres::AutoDiffCostFunction<ProjectErrorCostFunctionPinehole_SFM, 2, 3, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 4, 3>
			(new ProjectErrorCostFunctionPinehole_SFM(img_pixel)));
	}
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
	const Eigen::Vector2d m_img_pixel;
};
class ProjectErrorCostFunctionPinehole_SFM_first
{
public:
	ProjectErrorCostFunctionPinehole_SFM_first(const Eigen::Vector2d img_pixel) :
		m_img_pixel(img_pixel) {}

	template <typename T>
	bool operator()(const T* obj_3d, const T* const camK_F, const T* const camK_C, const T* const camK_shear,
		const T* const camDK_1, const T* const camDK_2, const T* const camDK_3, const T* const camDK_4, const T* const camDK_5, const T* const camDK_6,
		const T* const camDP, const T* const camDT_2, const T* const camDT_4, T* residuals) const
	{
		// 上面的参数类型为数组

		Eigen::Matrix<T, 2, 1> img = m_img_pixel.cast<T>();
		Eigen::Matrix<T, 3, 1> obj_cam_coor ;
		obj_cam_coor << obj_3d[0], obj_3d[1], obj_3d[2];
		T a = obj_cam_coor(0) / obj_cam_coor(2);
		T b = obj_cam_coor(1) / obj_cam_coor(2);

		T r2 = (a * a + b * b);
		T r4 = r2 * r2;
		T r6 = r2 * r4;

		T xd = a * (1.0 + camDK_1[0] * r2 + camDK_2[0] * r4 + camDK_3[0] * r6)
			/ (1.0 + camDK_4[0] * r2 + camDK_5[0] * r4 + camDK_6[0] * r6) +
			2.0 * camDP[0] * a * b + camDP[1] * (r2 + 2.0 * a * a) + camDT_2[0] * r2 + +camDT_4[0] * r4;
		T yd = b * (1.0 + camDK_1[0] * r2 + camDK_2[0] * r4 + camDK_3[0] * r6)
			/ (1.0 + camDK_4[0] * r2 + camDK_5[0] * r4 + camDK_6[0] * r6) +
			2.0 * camDP[1] * a * b + camDP[0] * (r2 + 2.0 * b * b) + camDT_2[1] * r2 + +camDT_4[1] * r4;

		T ud = camK_F[0] * xd + camK_shear[0] * yd + camK_C[0];
		T vd = camK_F[1] * yd + camK_C[1];

		// 残差 二维
		residuals[0] = ud - img(0);
		residuals[1] = vd - img(1);

		return true;
	}

	// 生成误差函数
	static ceres::CostFunction* Create(Eigen::Vector2d img_pixel)
	{
		return (new ceres::AutoDiffCostFunction<ProjectErrorCostFunctionPinehole_SFM_first, 2, 3, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2>
			(new ProjectErrorCostFunctionPinehole_SFM_first(img_pixel)));
	}
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
	const Eigen::Vector2d m_img_pixel;
};
class ProjectErrorCostFunctionPinehole_SFM_Rtheta
{
public:
	ProjectErrorCostFunctionPinehole_SFM_Rtheta(const Eigen::Vector2d img_pixel) :
		m_img_pixel(img_pixel) {}

	template <typename T>
	bool operator()(const T* obj_3d, const T* const camK_F, const T* const camK_C, const T* const camK_shear,
		const T* const camDK_1, const T* const camDK_2, const T* const camDK_3, const T* const camDK_4, const T* const camDK_5, const T* const camDK_6,
		const T* const camDP, const T* const camDT_2, const T* const camDT_4,
		const T* const camQvec, const T* const camTvec, T* residuals) const
	{
		// 上面的参数类型为数组

		Eigen::Matrix<T, 3, 1> obj;
		obj << obj_3d[0], obj_3d[1], obj_3d[2];


		Eigen::Quaternion<T> qvec(camQvec);
		Eigen::Matrix<T, 3, 1> tvec;
		tvec << sin(camTvec[0]), cos(camTvec[0])* sin(camTvec[1]), cos(camTvec[0])* cos(camTvec[1]);

		Eigen::Matrix<T, 2, 1> img = m_img_pixel.cast<T>();

		Eigen::Matrix<T, 3, 1> obj_cam_coor = qvec.toRotationMatrix() * obj + tvec;

		T a = obj_cam_coor(0) / obj_cam_coor(2);
		T b = obj_cam_coor(1) / obj_cam_coor(2);

		T r2 = (a * a + b * b);
		T r4 = r2 * r2;
		T r6 = r2 * r4;

		T xd = a * (1.0 + camDK_1[0] * r2 + camDK_2[0] * r4 + camDK_3[0] * r6)
			/ (1.0 + camDK_4[0] * r2 + camDK_5[0] * r4 + camDK_6[0] * r6) +
			2.0 * camDP[0] * a * b + camDP[1] * (r2 + 2.0 * a * a) + camDT_2[0] * r2 + +camDT_4[0] * r4;
		T yd = b * (1.0 + camDK_1[0] * r2 + camDK_2[0] * r4 + camDK_3[0] * r6)
			/ (1.0 + camDK_4[0] * r2 + camDK_5[0] * r4 + camDK_6[0] * r6) +
			2.0 * camDP[1] * a * b + camDP[0] * (r2 + 2.0 * b * b) + camDT_2[1] * r2 + +camDT_4[1] * r4;

		T ud = camK_F[0] * xd + camK_shear[0] * yd + camK_C[0];
		T vd = camK_F[1] * yd + camK_C[1];

		// 残差 二维
		residuals[0] = ud - img(0);
		residuals[1] = vd - img(1);

		return true;
	}

	// 生成误差函数
	static ceres::CostFunction* Create(Eigen::Vector2d img_pixel)
	{
		return (new ceres::AutoDiffCostFunction<ProjectErrorCostFunctionPinehole_SFM_Rtheta, 2, 3, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 4, 2>
			(new ProjectErrorCostFunctionPinehole_SFM_Rtheta(img_pixel)));
	}
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
	const Eigen::Vector2d m_img_pixel;
};


class ProjectErrorCostFunctionPinehole_SFM_sameF
{
public:
	ProjectErrorCostFunctionPinehole_SFM_sameF(const Eigen::Vector2d img_pixel) :
		m_img_pixel(img_pixel) {}

	template <typename T>
	bool operator()(const T* obj_3d, const T* const camK_F, const T* const camK_C, const T* const camK_shear,
		const T* const camDK_1, const T* const camDK_2, const T* const camDK_3,
		const T* const camDK_4, const T* const camDK_5, const T* const camDK_6,
		const T* const camDP, const T* const camDT_2, const T* const camDT_4,
		const T* const camQvec, const T* const camTvec, T* residuals) const
	{
		// 上面的参数类型为数组

		Eigen::Matrix<T, 3, 1> obj;
		obj << obj_3d[0], obj_3d[1], obj_3d[2];

		Eigen::Quaternion<T> qvec(camQvec);
		Eigen::Matrix<T, 3, 1> tvec;
		tvec << camTvec[0], camTvec[1], camTvec[2];

		Eigen::Matrix<T, 2, 1> img = m_img_pixel.cast<T>();

		Eigen::Matrix<T, 3, 1> obj_cam_coor = qvec.toRotationMatrix() * obj + tvec;

		T a = obj_cam_coor(0) / obj_cam_coor(2);
		T b = obj_cam_coor(1) / obj_cam_coor(2);

		T r2 = (a * a + b * b);
		T r4 = r2 * r2;
		T r6 = r2 * r4;

		T xd = a * (1.0 + camDK_1[0] * r2 + camDK_2[0] * r4 + camDK_3[0] * r6)
			/ (1.0 + camDK_4[0] * r2 + camDK_5[0] * r4 + camDK_6[0] * r6) +
			2.0 * camDP[0] * a * b + camDP[1] * (r2 + 2.0 * a * a) + camDT_2[0] * r2 + +camDT_4[0] * r4;
		T yd = b * (1.0 + camDK_1[0] * r2 + camDK_2[0] * r4 + camDK_3[0] * r6)
			/ (1.0 + camDK_4[0] * r2 + camDK_5[0] * r4 + camDK_6[0] * r6) +
			2.0 * camDP[1] * a * b + camDP[0] * (r2 + 2.0 * b * b) + camDT_2[1] * r2 + +camDT_4[1] * r4;

		T ud = camK_F[0] * xd + camK_shear[0] * yd + camK_C[0];
		T vd = camK_F[0] * yd + camK_C[1];

		// 残差 二维
		residuals[0] = ud - img(0);
		residuals[1] = vd - img(1);

		return true;
	}

	// 生成误差函数
	static ceres::CostFunction* Create(Eigen::Vector2d img_pixel)
	{
		return (new ceres::AutoDiffCostFunction<ProjectErrorCostFunctionPinehole_SFM_sameF, 2, 3, 1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 4, 3>
			(new ProjectErrorCostFunctionPinehole_SFM_sameF(img_pixel)));
	}
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
	const Eigen::Vector2d m_img_pixel;
};
class ProjectErrorCostFunctionPinehole_SFM_first_sameF
{
public:
	ProjectErrorCostFunctionPinehole_SFM_first_sameF(const Eigen::Vector2d img_pixel) :
		m_img_pixel(img_pixel) {}

	template <typename T>
	bool operator()(const T* obj_3d, const T* const camK_F, const T* const camK_C, const T* const camK_shear,
		const T* const camDK_1, const T* const camDK_2, const T* const camDK_3, const T* const camDK_4, const T* const camDK_5, const T* const camDK_6,
		const T* const camDP, const T* const camDT_2, const T* const camDT_4, T* residuals) const
	{
		// 上面的参数类型为数组

		Eigen::Matrix<T, 2, 1> img = m_img_pixel.cast<T>();
		Eigen::Matrix<T, 3, 1> obj_cam_coor;
		obj_cam_coor << obj_3d[0], obj_3d[1], obj_3d[2];
		T a = obj_cam_coor(0) / obj_cam_coor(2);
		T b = obj_cam_coor(1) / obj_cam_coor(2);

		T r2 = (a * a + b * b);
		T r4 = r2 * r2;
		T r6 = r2 * r4;

		T xd = a * (1.0 + camDK_1[0] * r2 + camDK_2[0] * r4 + camDK_3[0] * r6)
			/ (1.0 + camDK_4[0] * r2 + camDK_5[0] * r4 + camDK_6[0] * r6) +
			2.0 * camDP[0] * a * b + camDP[1] * (r2 + 2.0 * a * a) + camDT_2[0] * r2 + +camDT_4[0] * r4;
		T yd = b * (1.0 + camDK_1[0] * r2 + camDK_2[0] * r4 + camDK_3[0] * r6)
			/ (1.0 + camDK_4[0] * r2 + camDK_5[0] * r4 + camDK_6[0] * r6) +
			2.0 * camDP[1] * a * b + camDP[0] * (r2 + 2.0 * b * b) + camDT_2[1] * r2 + +camDT_4[1] * r4;

		T ud = camK_F[0] * xd + camK_shear[0] * yd + camK_C[0];
		T vd = camK_F[0] * yd + camK_C[1];

		// 残差 二维
		residuals[0] = ud - img(0);
		residuals[1] = vd - img(1);

		return true;
	}

	// 生成误差函数
	static ceres::CostFunction* Create(Eigen::Vector2d img_pixel)
	{
		return (new ceres::AutoDiffCostFunction<ProjectErrorCostFunctionPinehole_SFM_first_sameF, 2, 3, 1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2>
			(new ProjectErrorCostFunctionPinehole_SFM_first_sameF(img_pixel)));
	}
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
	const Eigen::Vector2d m_img_pixel;
};
class ProjectErrorCostFunctionPinehole_SFM_Rtheta_sameF
{
public:
	ProjectErrorCostFunctionPinehole_SFM_Rtheta_sameF(const Eigen::Vector2d img_pixel) :
		m_img_pixel(img_pixel) {}

	template <typename T>
	bool operator()(const T* obj_3d, const T* const camK_F, const T* const camK_C, const T* const camK_shear,
		const T* const camDK_1, const T* const camDK_2, const T* const camDK_3, const T* const camDK_4, const T* const camDK_5, const T* const camDK_6,
		const T* const camDP, const T* const camDT_2, const T* const camDT_4,
		const T* const camQvec, const T* const camTvec, T* residuals) const
	{
		// 上面的参数类型为数组

		Eigen::Matrix<T, 3, 1> obj;
		obj << obj_3d[0], obj_3d[1], obj_3d[2];


		Eigen::Quaternion<T> qvec(camQvec);
		Eigen::Matrix<T, 3, 1> tvec;
		tvec << sin(camTvec[0]), cos(camTvec[0])* sin(camTvec[1]), cos(camTvec[0])* cos(camTvec[1]);

		Eigen::Matrix<T, 2, 1> img = m_img_pixel.cast<T>();

		Eigen::Matrix<T, 3, 1> obj_cam_coor = qvec.toRotationMatrix() * obj + tvec;

		T a = obj_cam_coor(0) / obj_cam_coor(2);
		T b = obj_cam_coor(1) / obj_cam_coor(2);

		T r2 = (a * a + b * b);
		T r4 = r2 * r2;
		T r6 = r2 * r4;

		T xd = a * (1.0 + camDK_1[0] * r2 + camDK_2[0] * r4 + camDK_3[0] * r6)
			/ (1.0 + camDK_4[0] * r2 + camDK_5[0] * r4 + camDK_6[0] * r6) +
			2.0 * camDP[0] * a * b + camDP[1] * (r2 + 2.0 * a * a) + camDT_2[0] * r2 + +camDT_4[0] * r4;
		T yd = b * (1.0 + camDK_1[0] * r2 + camDK_2[0] * r4 + camDK_3[0] * r6)
			/ (1.0 + camDK_4[0] * r2 + camDK_5[0] * r4 + camDK_6[0] * r6) +
			2.0 * camDP[1] * a * b + camDP[0] * (r2 + 2.0 * b * b) + camDT_2[1] * r2 + +camDT_4[1] * r4;

		T ud = camK_F[0] * xd + camK_shear[0] * yd + camK_C[0];
		T vd = camK_F[0] * yd + camK_C[1];

		// 残差 二维
		residuals[0] = ud - img(0);
		residuals[1] = vd - img(1);

		return true;
	}

	// 生成误差函数
	static ceres::CostFunction* Create(Eigen::Vector2d img_pixel)
	{
		return (new ceres::AutoDiffCostFunction<ProjectErrorCostFunctionPinehole_SFM_Rtheta_sameF, 2, 3, 1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 4, 2>
			(new ProjectErrorCostFunctionPinehole_SFM_Rtheta_sameF(img_pixel)));
	}
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
	const Eigen::Vector2d m_img_pixel;
};

class SFM_optimize
{
public:
	static void calculate_SFM(std::vector<std::vector<Coded_detect_inf>> code_circle_serial, 
		std::vector<sfm_3d_group> &point_clouds,
		std::vector<Eigen::Quaterniond>& camQvec, std::vector<Eigen::Vector3d>& camTvec
		, std::vector<std::vector<cv::Point2d>>& Re_project_Map
		, double*& camK, double*& camDK, double*& camDP, double*& camDT
		, ceres::Solver::Summary* summary
		, ceres::Solver::Summary* summary_AG
		, unsigned int max_iter_num = 1000
		, double stop_value = 1e-9
		, unsigned int num_threads = 1
		, unsigned int timeout = 1e8
		, unsigned char Loss_type = 3
		, double Loss_value = 1
		, unsigned int Dis_K_num = 3
		, unsigned int Dis_P_num = 2
		, unsigned int Dis_T_num = 0
		, bool Fixed_all = false
		, bool Use_F = true
		, bool Use_Cx_Cy = true
		, bool Use_shear = false
		, bool Use_same_F = false
		, double* F_range = nullptr
		, double* C_range = nullptr
		, double out_th = 10);
	static void calculate_dual(std::vector<std::vector<Coded_detect_inf>> code_circle_serial,
		std::vector<sfm_3d_group>& point_clouds,
		std::vector<Eigen::Quaterniond>& camQvec, std::vector<Eigen::Vector3d>& camTvec
		, std::vector<std::vector<cv::Point2d>>& Re_project_Map
		, double*& camK, double*& camDK, double*& camDP, double*& camDT
		, ceres::Solver::Summary* summary
		, unsigned int max_iter_num = 1000
		, double stop_value = 1e-9
		, unsigned int num_threads = 1
		, unsigned int timeout = 1e8
		, unsigned char Loss_type = 3
		, double Loss_value = 1
		, unsigned int Dis_K_num = 3
		, unsigned int Dis_P_num = 2
		, unsigned int Dis_T_num = 0
		, bool Fixed_all = false
		, bool Use_F = true
		, bool Use_Cx_Cy = true
		, bool Use_shear = false
		, bool Use_same_F = false
		, double* F_range = nullptr
		, double* C_range = nullptr);

	static void calculate_RT(std::vector<cv::Point3d> point_clouds,
		std::vector<cv::Point2d> img_points,
		Eigen::Quaterniond& camQvec, Eigen::Vector3d& camTvec
		, double* camK, double* camDK, double* camDP, double* camDT
		, unsigned int max_iter_num = 1000
		, double stop_value = 1e-9
		, unsigned int num_threads = 1
		, unsigned int timeout = 1e8
		, unsigned char Loss_type = 3
		, double Loss_value = 1);

	static void calculate_Circle(std::vector<sfm_3d_group> &point_clouds, 
		std::vector<std::vector<std::vector<cv::Point2f>>> contour_circle_serial,
		std::vector<std::vector<Coded_detect_inf>> code_circle_serial, std::vector<Eigen::Vector3d> &ori_serial,
		std::vector<Eigen::Quaterniond> camQvec, std::vector<Eigen::Vector3d> camTvec
		, double* camK, double* camDK, double* camDP, double* camDT
		, std::vector<ceres::Solver::Summary>* summary
		, unsigned int max_iter_num = 1000
		, double stop_value = 1e-9
		, unsigned int num_threads = 1
		, unsigned int timeout = 1e8
		, unsigned char Loss_type = 3
		, double Loss_value = 1);
};

