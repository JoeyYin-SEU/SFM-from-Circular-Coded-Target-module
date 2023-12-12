#pragma once
/*支持的相机类型*/
#include <QString>
#define Max_sectorSize 4096
#define Max_Camera_support_Disk_size 4
#define maxFileNr 10

#define Max_camera_support 25
#define Max_bining_skipping_size 64
#define Max_shutter_mode_size 16
#define Max_quality_size 1000
#define Max_Disk_size 6

/*自定义相机参数存储结构*/
struct Camera_inf
{
	Camera_inf() :Cam_serial(""), frame_rate(0), frame_rate_max(0), frame_rate_min(0), exposure_time(0), exposure_time_max(0),
		exposure_time_min(0), gain(0), gain_max(0), gain_min(0), gamma(0), gamma_max(0), gamma_min(0),
		width_pixel(0), width_pixel_max(0), width_pixel_min(0), height_pixel(0), height_pixel_max(0),
		height_pixel_min(0), offx_pixel(-1), offx_pixel_max(-1), offx_pixel_min(-1), offy_pixel(-1), offy_pixel_max(-1),
		offy_pixel_min(-1), tempera_1(0), tempera_2(0), width_pixel_inc(0), height_pixel_inc(0), offx_pixel_inc(-1),
		offy_pixel_inc(-1), bin_V_now(-1), bin_H_now(-1), skip_V_now(-1), skip_H_now(-1), trigger_mode(-2), shutter_mode_now(-1), skip_stream_number(1)
	{
		for (int ii = 0; ii < Max_bining_skipping_size; ii++)
		{
			bin_V[ii] = -1;
			bin_H[ii] = -1;
			skip_V[ii] = -1;
			skip_H[ii] = -1;
		}
		for (int ii = 0; ii < Max_shutter_mode_size; ii++)
		{
			shutter_mode[ii] = -1;
		}
	}
	QString Cam_serial;
	double frame_rate;
	double frame_rate_max;
	double frame_rate_min;
	double exposure_time;
	double exposure_time_max;
	double exposure_time_min;
	double gain;
	double gain_max;
	double gain_min;
	double gamma;
	double gamma_max;
	double gamma_min;
	int width_pixel;
	int width_pixel_max;
	int width_pixel_min;
	int width_pixel_inc;
	int height_pixel;
	int height_pixel_max;
	int height_pixel_min;
	int height_pixel_inc;
	int offx_pixel;
	int offx_pixel_max;
	int offx_pixel_min;
	int offx_pixel_inc;
	int offy_pixel;
	int offy_pixel_max;
	int offy_pixel_min;
	int offy_pixel_inc;
	float tempera_1;
	float tempera_2;
	int bin_V[Max_bining_skipping_size];
	int bin_H[Max_bining_skipping_size];
	int bin_V_now;
	int bin_H_now;
	int skip_V[Max_bining_skipping_size];
	int skip_H[Max_bining_skipping_size];
	int skip_V_now;
	int skip_H_now;
	int trigger_mode;
	int shutter_mode[Max_shutter_mode_size];
	int shutter_mode_now;
	int skip_stream_number;
};

struct DIC_Min_Max_params
{
	double max_value;
	double min_value;
	DIC_Min_Max_params()
	{
	}
};

enum DIC_auto_ROI_type
{
	Grad_std,
	ZNCC,
};

struct DIC_auto_ROI_params
{
	int window_R;
	double threshold;
	double threshold_2;
	int Step;
	DIC_auto_ROI_params()
	{
	}
};
enum DIC_failed_type
{
	Out_of_borded,
	Out_of_init_zncc,
	Out_of_icgn_error,
	No_texture,
	Gradient_explosion
};
enum DIC_order
{
	First_order,
	Second_order
};

enum Window_Cross_type
{
	Circle_type,
	Square_type
};
enum Interpolation_method
{
	Neareast_Interpolation,
	Liner_Interpolation,
	Cubic_Interpolation,
	Quintic_Interpolation
};
struct ellipse_view_pose
{
	double view_N_1[3];
	double view_N_2[3];
	ellipse_view_pose()
	{
		view_N_1[0] = 0;
		view_N_1[1] = 0;
		view_N_1[2] = 0;
		view_N_2[0] = 0;
		view_N_2[1] = 0;
		view_N_2[2] = 0;
	}
};

struct DIC_Params_double
{
	int radius;
	int interval;
	bool used_gpu;
	Window_Cross_type area_type;
	Interpolation_method interpolation_base;
	DIC_order order;
	int search_R;
	float Min_corr;
	float Min_error;
	float stop_g;
	int max_iter;
	DIC_Params_double()
	{
		radius = 11;
		interval = 1;
		used_gpu = false;
		area_type = Window_Cross_type::Circle_type;
		interpolation_base = Interpolation_method::Quintic_Interpolation;
		order = DIC_order::First_order;
		search_R = 10;
		Min_corr = 0.7;
		Min_error = 1e-2;
		stop_g = 1e-4;
		max_iter = 10;
	}
};
struct DIC_Params
{
	int radius;
	int interval;
	bool used_gpu;
	Window_Cross_type area_type;
	Interpolation_method interpolation_base;
	DIC_order order;
	int search_R;
	float Min_corr;
	float Min_error; 
	float stop_g; 
	int max_iter;
	DIC_Params()
	{
		radius = 11;
		interval = 1;
		used_gpu = false;
		area_type = Window_Cross_type::Circle_type;
		interpolation_base = Interpolation_method::Quintic_Interpolation;
		order = DIC_order::First_order;
		search_R = 10;
		Min_corr = 0.7;
		Min_error = 1e-2;
		stop_g = 1e-4;
		max_iter = 10;
	}
};

struct match_3d_point_group
{
	int matching_number;
	int position;
	match_3d_point_group(int val, int x) :matching_number(val), position(x)
	{
	}
	/*匹配点越高越好*/
	friend bool operator < (struct match_3d_point_group a, struct match_3d_point_group b)
	{
		return a.matching_number < b.matching_number;
	}
};

struct match_pair_group
{
	int matching_number;
	int position_x;
	int position_y;
	match_pair_group(int val, int x, int y) :matching_number(val), position_x(x), position_y(y)
	{
	}
	/*匹配点越高越好*/
	friend bool operator < (struct match_pair_group a, struct match_pair_group b)
	{
		return a.matching_number < b.matching_number;
	}
};
struct sfm_3d_group
{
	int code;
	int weight;
	double x;
	double y;
	double z;
	double R;
	double error;
	sfm_3d_group(int val, int w, double x_val, double y_val, double z_val, double err, double RS = 0) :code(val), weight(w), x(x_val), y(y_val), z(z_val), error(err), R(RS)
	{

	}
	/*误差越低越好*/
	friend bool operator < (struct sfm_3d_group a, struct sfm_3d_group b)
	{
		return a.error >= b.error;
	}
};

struct Coded_detect_inf
{
	Coded_detect_inf() :x(0), y(0), r_a(0), r_b(0), angle_in_pi(0), code_num(-1), fit_err(-1), image_index(-1)
	{
	}
	double x;
	double y;
	double r_a;
	double r_b;
	double angle_in_pi;
	int code_num;
	double fit_err;
	int image_index;
	ellipse_view_pose pose;
};

enum Cameras_type
{
	Hik_camera_type,
	Flir_camera_type,
	Ids_camera_type,
	Kaya_camera_type
};
/*支持的压电陶瓷类型*/
enum PZT_type
{
	Thorlabs_KPZ101_type
};
/*支持的投影仪类型*/
enum Projector_type
{
	DLP6500_type,
	DLP4500_type
};
/*支持的触发器类型*/
enum TriggerControl_type
{
	GuanLi_type
};
/*支持的平移旋转台类型*/
enum Tran_Roll_type
{
	SHOT702_80_20ZF_type
};

enum Target_type
{
	Chess_Board_type,
	Circle_Board_type,
	Ori_Circle_Board_type,
	Speckle_Board_type
};

/*数据类型
* 0:图像数据
* 1:标定角点数据
* 2:标定参数数据
* 3:点云数据
* other:无效数据
*/

/*数据格式
* 0:unsigned char 1字节
* 1:float 4字节
* 2:double 8字节
* other:无效数据
*/
/*MSDATA格式组成(图像格式下)：
*unsigned char-数据类型 1字节 = 0
*char-数据格式 1字节
*bool-存在软索引 1字节
*bool-存在硬索引 1字节
*bool-存在软时间 1字节
*bool-存在硬时间 1字节
*bool-存在温度-1 1字节
*bool-存在温度-2 1字节
*bool-备用-1 1字节
*bool-备用-2 1字节
*bool-备用-3 1字节
*bool-备用-4 1字节
*bool-备用-5 1字节
*bool-备用-6 1字节
*bool-备用-7 1字节
*bool-备用-8 1字节
*bool-备用-9 1字节
*bool-存在图像 1字节
* 补齐对齐字节
** for以下循环写入
* 图像宽度 int 4字节
* 图像高度 int 4字节
* 图像软索引（如果索引TRUE） long long 8字节
* 图像硬索引（如果索引TRUE） long long 8字节
* 图像软时间（如果索引TRUE） long long 8字节
* 图像硬时间（如果索引TRUE） long long 8字节
* 图像温度-1（如果索引TRUE） float 4字节
* 图像温度-2（如果索引TRUE） float 4字节
* 备用-1
* 备用-2
* 备用-3
* 备用-4
* 备用-5
* 备用-6
* 备用-7
* 备用-8
* 图像数据（如果索引TRUE） 图像高度*图像深度*图像位数字节
* 补齐对齐字节
** end
*/

/*MSDATA格式组成(Chess标定角点格式下)：
*unsigned char-数据类型 1字节 = 1(chess:1/circle:3)
** for以下循环写入
* 相机组号 int 4字节
* 相机序号 int 4字节
* 图像序号 int 4字节
* 标定图像位置字节大小 long long 8字节
* 标定图像位置 char 由字节大小读取
* 图像宽度 int 4字节
* 图像高度 int 4字节
* 标定角点数量 int 4字节
* 角点行数 int 4字节
* 角点列数 int 4字节
**** for以下循环写入角点
* X轴像素位置 double 8字节
* Y轴像素位置 double 8字节
**** end
** end
*/

// 标定结果格式的数据内容包含角点，属于更高层次的保存
/*MSDATA格式组成(标定结果格式下)：
*unsigned char-数据类型 1字节 = 2(chess:2/circle:4)
** 以下循环写入
* 相机组号 int 4字节
* 相机序号 int 4字节
* 相机内参 double 8字节*5 [fx,fy,cx,cy,s] (或等价于此排布的随机矩阵形式)
* 相机内参K畸变 double 8字节*6 [K1,K2,K3,K4,K5,K6]
* 相机内参P畸变 double 8字节*2 [P1,P2]
* 相机内参Thin畸变 double 8字节*4 [T1,T2,T3,T4]
* 相机外参旋转四元数 double 8字节*4 [x,y,z,w] (或等价于此排布的随机矩阵形式)
* 相机外参平移矩阵 double 8字节*3 [Tx,Ty,Tz] (或等价于此排布的随机矩阵形式)
* 验证组数 int 4字节
*** for以下循环写入（循环[验证组数]次）
* 相机内参 double 8字节*5 [fx,fy,cx,cy,s] (或等价于此排布的随机矩阵形式)
* 相机内参K畸变 double 8字节*6 [K1,K2,K3,K4,K5,K6]
* 相机内参P畸变 double 8字节*2 [P1,P2]
* 相机内参Thin畸变 double 8字节*4 [T1,T2,T3,T4]
* 相机外参旋转四元数 double 8字节*4 [x,y,z,w] (或等价于此排布的随机矩阵形式)
* 相机外参平移矩阵 double 8字节*3 [Tx,Ty,Tz] (或等价于此排布的随机矩阵形式)
* 角点数量 int 4字节
* (矫正后)标定板形貌 [x,y,z] double 8字节*3*角点数量
*** end
* 相机重投影 double 8字节
* 单相机迭代次数 int 4字节
* 单相机迭代代价 double 4字节*单相机迭代次数
* 多相机迭代次数 int 4字节
* 多相机迭代代价 double 4字节*多相机迭代次数
* 角点数量 int 4字节
* (矫正后)标定板形貌 [x,y,z] double 8字节*3*角点数量
* 相机图像数量 int 4字节
*** for以下循环写入
* 图像序号 int 4字节
* 图像使能 bool 1字节
* 图像重投影误差 double 8字节
* 单视图旋转四元数 double 8字节*4
* 单视图平移矩阵 double 8字节*3
* 标定图像位置字节大小 long long 8字节
* 标定图像位置 char 由字节大小读取
* 图像宽度 int 4字节
* 图像高度 int 4字节
* 标定角点数量 int 4字节
* 角点行数 int 4字节
* 角点列数 int 4字节
**** for以下循环写入角点信息
* X轴像素位置 double 8字节
* Y轴像素位置 double 8字节
* X轴像素误差 double 8字节
* Y轴像素误差 double 8字节
**** end
*** end
*/


/*MSDATA格式组成(oricirle标定角点格式下)：
*unsigned char-数据类型 1字节 = 5
** for以下循环写入
* 相机组号 int 4字节
* 相机序号 int 4字节
* 图像序号 int 4字节
* 标定图像位置字节大小 long long 8字节
* 标定图像位置 char 由字节大小读取
* 图像宽度 int 4字节
* 图像高度 int 4字节
* 标定角点数量 int 4字节
* 角点行数 int 4字节
* 角点列数 int 4字节
* 角点offset行数 int 4字节
* 角点offset列数 int 4字节
* 角点incre行数 int 4字节
* 角点incre列数 int 4字节
**** for以下循环写入角点
* X轴像素位置 double 8字节
* Y轴像素位置 double 8字节
**** end
** end
*/

// 标定结果格式的数据内容包含角点，属于更高层次的保存
/*MSDATA格式组成(oricirle标定结果格式下)：
*unsigned char-数据类型 1字节 = 6
** 以下循环写入
* 相机组号 int 4字节
* 相机序号 int 4字节
* 相机内参 double 8字节*5 [fx,fy,cx,cy,s] (或等价于此排布的随机矩阵形式)
* 相机内参K畸变 double 8字节*6 [K1,K2,K3,K4,K5,K6]
* 相机内参P畸变 double 8字节*2 [P1,P2]
* 相机内参Thin畸变 double 8字节*4 [T1,T2,T3,T4]
* 相机外参旋转四元数 double 8字节*4 [x,y,z,w] (或等价于此排布的随机矩阵形式)
* 相机外参平移矩阵 double 8字节*3 [Tx,Ty,Tz] (或等价于此排布的随机矩阵形式)
* 验证组数 int 4字节
*** for以下循环写入（循环[验证组数]次）
* 相机内参 double 8字节*5 [fx,fy,cx,cy,s] (或等价于此排布的随机矩阵形式)
* 相机内参K畸变 double 8字节*6 [K1,K2,K3,K4,K5,K6]
* 相机内参P畸变 double 8字节*2 [P1,P2]
* 相机内参Thin畸变 double 8字节*4 [T1,T2,T3,T4]
* 相机外参旋转四元数 double 8字节*4 [x,y,z,w] (或等价于此排布的随机矩阵形式)
* 相机外参平移矩阵 double 8字节*3 [Tx,Ty,Tz] (或等价于此排布的随机矩阵形式)
* 角点数量 int 4字节
* (矫正后)标定板形貌 [x,y,z] double 8字节*3*角点数量
*** end
* 相机重投影 double 8字节
* 单相机迭代次数 int 4字节
* 单相机迭代代价 double 4字节*单相机迭代次数
* 多相机迭代次数 int 4字节
* 多相机迭代代价 double 4字节*多相机迭代次数
* 角点数量 int 4字节
* (矫正后)标定板形貌 [x,y,z] double 8字节*3*角点数量
* 相机图像数量 int 4字节
*** for以下循环写入
* 图像序号 int 4字节
* 图像使能 bool 1字节
* 图像重投影误差 double 8字节
* 单视图旋转四元数 double 8字节*4
* 单视图平移矩阵 double 8字节*3
* 标定图像位置字节大小 long long 8字节
* 标定图像位置 char 由字节大小读取
* 图像宽度 int 4字节
* 图像高度 int 4字节
* 标定角点数量 int 4字节
* 角点行数 int 4字节
* 角点列数 int 4字节
* 角点offset行数 int 4字节
* 角点offset列数 int 4字节
* 角点incre行数 int 4字节
* 角点incre列数 int 4字节
**** for以下循环写入角点信息
* X轴像素位置 double 8字节
* Y轴像素位置 double 8字节
* X轴像素误差 double 8字节
* Y轴像素误差 double 8字节
**** end
*** end
*/

/*MSDATA格式组成(SPSI结果格式下)：
*unsigned char-数据类型 1字节 = 21
* 图像宽度 int 4字节
* 图像高度 int 4字节
**** for 以下循环写入相位
* 相位值 float 4字节
**** end
*/

/*MSDATA格式组成(摄影测量编码圆特征点格式下)：
*unsigned char-数据类型 1字节 = 31
** for以下循环写入
* 图像序号 int 4字节
* 图像位置字节大小 long long 8字节
* 图像位置 char 由字节大小读取
* 图像宽度 int 4字节
* 图像高度 int 4字节
* 特征点数量 int 4字节
**** for以下循环写入角点
* X轴像素位置 double 8字节
* Y轴像素位置 double 8字节
* 长轴长度 double 8字节
* 短轴长度 double 8字节
* 旋转角 double 8字节
* 编码值 int 4字节
* 拟合误差 double 8字节
* 圆点标志物法向猜测-1 double 8字节*3
* 圆点标志物法向猜测-2 double 8字节*3
* 轮廓点数 int 4字节
* 轮廓值[x,y] 2*int 4字节*轮廓点数
**** end
** end
*/


/*MSDATA格式组成(摄影测量编码圆结果格式下)：
* unsigned char-数据类型 1字节 = 32
* 相机内参 double 8字节*5 [fx,fy,cx,cy,s] (或等价于此排布的随机矩阵形式)
* 相机内参K畸变 double 8字节*6 [K1,K2,K3,K4,K5,K6]
* 相机内参P畸变 double 8字节*2 [P1,P2]
* 相机内参Thin畸变 double 8字节*4 [T1,T2,T3,T4]
* 图像数量 int 4字节
** for以下循环写入
* 图像序号 int 4字节
* 图像位置字节大小 long long 8字节
* 图像位置 char 由字节大小读取
* 图像宽度 int 4字节
* 图像高度 int 4字节
* 单视图旋转四元数 8字节*4
* 单视图平移矩阵 8字节*3
* 特征点数量 int 4字节
**** for以下循环写入角点
* X轴像素位置 double 8字节
* Y轴像素位置 double 8字节
* 长轴长度 double 8字节
* 短轴长度 double 8字节
* 旋转角 double 8字节
* 编码值 int 4字节
* 拟合误差 double 8字节
* 圆点标志物法向猜测-1 double 8字节*3
* 圆点标志物法向猜测-2 double 8字节*3
* 轮廓点数 int 4字节
* 轮廓值[x,y] 2*int 4字节*轮廓点数
* 重投影误差 double 8字节*2
**** end
** end
* 三维点数量 int 4字节
** for以下循环写入三维点
* 编码值 int 4字节
* 视图权重 int 4字节
* X轴像素位置 double 8字节
* Y轴像素位置 double 8字节
* Z轴像素位置 double 8字节
* 圆点标志物法向猜测 double 8字节*3
** end
* 计算序列数量 int 4字节
** for以下循环写入计算序列
* 序列索引 int 4字节
** end
* 尺度信息数量 int 4字节
** for以下循环写入尺度
* 起始编码点编码值 int 4字节
* 终止编码点编码值 int 4字节
* 距离值 double 8 字节
** end
* 校平信息数量 int 4字节
** for以下循环写入尺度
* 需校平编码点编码值 int 4字节
** end
*/