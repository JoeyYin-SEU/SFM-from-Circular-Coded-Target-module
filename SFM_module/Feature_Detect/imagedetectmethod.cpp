#include "imagedetectmethod.h"
#include <imgproc/types_c.h>

ImageDetectMethod::ImageDetectMethod(QObject* parent)
	: QObject(parent)
{
}

ImageDetectMethod::~ImageDetectMethod()
{
}

//////**********************************2014-2-9
bool ImageDetectMethod::CodeAndUncodePointDetect(QString image_file_name, cv::Mat& code_point_mat, cv::Mat& uncode_point_mat,
	float ratio_k, float ratio_k1, float ratio_k2,
	float min_radius, float max_radius, float ellipse_error_pixel /*=0.5*/,
	MarkPointColorType color_type /*= BlackDownWhiteUp*/, CodePointBitesType code_bites_type /*=CodeBites15*/,
	DetectContoursMethod image_process_method /*= CANNY_Method*/, SubPixelPosMethod subpixel_pos_method /*=Gray_Centroid*/)
{
	cv::Mat ori_image = cv::imread(image_file_name.toStdString(), 0);
	if (!ori_image.data)        // 判断图片调入是否成功
		return false;        // 调入图片失败则退出

	cv::Mat processed_image_mat;
	std::vector<std::vector<cv::Point>> contours;
	QList<QList<float>> ellipse_pars;      //ellipse_pars - n*6  center_x,center_y,r_a,r_b,angle_inPI,ellipse_error,contours_index,ID,code_type(0- uncode,1- code)

	//图像预处理
	ImagePreprocess(ori_image, processed_image_mat);//高斯滤波前处理，C

	//边缘检测，存入闭合轮廓
	DetectClosedContours(processed_image_mat, contours, image_process_method);

	//轮廓筛选，尺寸，形状等准则
	FilterEllipseContours(contours, min_radius, max_radius,
		ellipse_error_pixel, ellipse_pars);

	//进一步筛选 ，用于编码点和非编码
	FilterEllipseContoursForCodePoint(processed_image_mat, ratio_k, ratio_k1, ratio_k2,
		ellipse_pars);

	//test
	//ShowContous( ori_image.cols,ori_image.rows,contours,ellipse_pars );
	//

	//////编码点与非编码点进行区分,解码
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
			ratio_k1, ratio_k2, color_type, code_bites_type);

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
			//bool is_uncodepoint = UncodePointCheck(processed_image_mat, ellipse_pars[i][0], ellipse_pars[i][1], ellipse_pars[i][2], ellipse_pars[i][3], ellipse_pars[i][4],
				//ratio_k, color_type, code_bites_type);
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

	//************************亚像素定位
	//ellipse_pars - n*6  center_x,center_y,r_a,r_b,angle_inPI,ellipse_error,contours_index,ID,code_type(0- uncode,1- code)
	std::vector<std::vector<cv::Point2f>> subpixel_edge_contours;
	for (int i = 0; i < ellipse_pars.size(); i++)
	{
		float sub_pixel_x, sub_pixel_y;
		std::vector<cv::Point2f> edge_contour;
		FindSubPixelPosOfCircleCenter20140210(processed_image_mat, ellipse_pars[i][0], ellipse_pars[i][1], ellipse_pars[i][2],
			ellipse_pars[i][3], ellipse_pars[i][4], contours[ellipse_pars[i][6]], sub_pixel_x, sub_pixel_y,
			&edge_contour, subpixel_pos_method, color_type);
		ellipse_pars[i][0] = sub_pixel_x;
		ellipse_pars[i][1] = sub_pixel_y;
		subpixel_edge_contours.push_back(edge_contour);
	}

	//输出  code_point_mat,  uncode_point_mat
	code_point_mat = cv::Mat();
	uncode_point_mat = cv::Mat();
	for (int i = 0; i < ellipse_pars.size(); i++)
	{
		//id,x,y,quality, r_a,r_b,angle_in_pi
		float a[7] = { ellipse_pars[i][7],ellipse_pars[i][0],ellipse_pars[i][1],ellipse_pars[i][5],
			ellipse_pars[i][2],ellipse_pars[i][3],ellipse_pars[i][4] };
		cv::Mat mat = cv::Mat(1, 7, CV_32F, a);

		//*******测试用
// 		double a[3] ={ellipse_pars[i][0],ellipse_pars[i][1],ellipse_pars[i][7]};
// 		Mat mat = Mat(1,3,CV_64F,a);
		///

		if (ellipse_pars[i][8] > 0)
		{
			code_point_mat.push_back(mat);
		}
		else
		{
			uncode_point_mat.push_back(mat);
		}
	}

	for (int ii = 0; ii < code_point_mat.rows; ii++)
	{
		//std::cout << uncode_point_mat.at<float>(ii, 0) << endl;
		circle(ori_image, cv::Point(code_point_mat.at<float>(ii,1), code_point_mat.at<float>(ii,2)),5, cv::Scalar(255,0,0));
	}

	//cv::namedWindow("img", cv::WINDOW_NORMAL);
	//imshow("img", ori_image);
	//cv::waitKey(0);
	return true;
}
void ImageDetectMethod::_sobel_gradient(const cv::Mat& mat, cv::Mat& dx, cv::Mat& dy, cv::Mat& magnitudes, cv::Mat& angles,
	int apertureSize, bool L2gradient)
{
	CV_Assert(apertureSize == 3 || apertureSize == 5);

	double scale = 1.0;
	cv::Sobel(mat, dx, CV_16S, 1, 0, apertureSize, scale, cv::BORDER_REPLICATE);
	cv::Sobel(mat, dy, CV_16S, 0, 1, apertureSize, scale, cv::BORDER_REPLICATE);

	const int TAN225 = 13573;			//tan22.5 * 2^15(2 << 15)

	angles = cv::Mat(mat.size(), CV_8UC1);  // 0-> horizontal, 1 -> vertical, 2 -> diagonal
	magnitudes = cv::Mat::zeros(mat.rows + 2, mat.cols + 2, CV_32SC1);
	cv::Mat magROI = cv::Mat(magnitudes, cv::Rect(1, 1, mat.cols, mat.rows));

	for (int i = 0; i < mat.rows; i++)
	{
		for (int j = 0; j < mat.cols; j++)
		{
			short xs = dx.ptr<short>(i)[j];
			short ys = dy.ptr<short>(i)[j];
			int x = (int)std::abs(xs);
			int y = (int)std::abs(ys) << 15;

			if (L2gradient) {
				//magROI.ptr<int>(i)[j] = int(xs) * xs + int(ys) * ys;
				magROI.ptr<int>(i)[j] = (int)std::sqrt(xs * xs + ys * ys);
			}
			else {
				magROI.ptr<int>(i)[j] = std::abs(int(xs)) + std::abs(int(ys));
			}

			int tan225x = x * TAN225;
			if (y < tan225x) {  // horizontal
				angles.ptr<uchar>(i)[j] = 0;
			}
			else
			{
				int tan675x = tan225x + (x << 16);
				if (y > tan675x) {  // vertical
					angles.ptr<uchar>(i)[j] = 1;
				}
				else {  // diagonal
					angles.ptr<uchar>(i)[j] = 2;
				}
			}
		}
	}
}
void ImageDetectMethod::_calculate_hysteresis_threshold_value(const cv::Mat& dx, const cv::Mat& dy, const cv::Mat& magnitudes,
	const cv::Mat& angles, cv::Mat& NMSImage, int& low, int& high)
{
	NMSImage = cv::Mat::zeros(magnitudes.size(), magnitudes.type());		//CV_32SC1

	for (int i = 0; i < dx.rows; ++i)
	{
		int r = i + 1;
		for (int j = 0; j < dx.cols; ++j)
		{
			int c = j + 1;
			int m = magnitudes.ptr<int>(r)[c];
			uchar angle = angles.ptr<uchar>(i)[j];

			if (angle == 0)			//horizontal
			{
				if (m > magnitudes.ptr<int>(r)[c - 1] && m >= magnitudes.ptr<int>(r)[c + 1])
					NMSImage.ptr<int>(r)[c] = m;
			}
			else if (angle == 1)	//vertical
			{
				if (m > magnitudes.ptr<int>(r - 1)[c] && m >= magnitudes.ptr<int>(r + 1)[c])
					NMSImage.ptr<int>(r)[c] = m;
			}
			else if (angle == 2)	//diagonal
			{
				short xs = dx.ptr<short>(i)[j];
				short ys = dy.ptr<short>(i)[j];
				if ((xs > 0 && ys > 0) || (xs < 0 && ys < 0))
				{	//45 degree
					if (m > magnitudes.ptr<int>(r - 1)[c - 1] && m > magnitudes.ptr<int>(r + 1)[c + 1])
						NMSImage.ptr<int>(r)[c] = m;
				}
				else
				{	//135 degree
					if (m > magnitudes.ptr<int>(r - 1)[c + 1] && m > magnitudes.ptr<int>(r + 1)[c - 1])
						NMSImage.ptr<int>(r)[c] = m;
				}
			}
		}
	}

	//利用Otsu对非极大值抑制图像进行处理，将计算得到的阈值作为高阈值high, 低阈值取高阈值的0.5倍
	cv::normalize(NMSImage, NMSImage, 0, 255, cv::NORM_MINMAX);
	NMSImage.convertTo(NMSImage, CV_8UC1);

	cv::Mat temp;
	high = (int)cv::threshold(NMSImage, temp, 0, 255, cv::THRESH_OTSU);
	low = (int)(0.5 * high);
}
void ImageDetectMethod::_non_maximum_suppression(const cv::Mat& NMSImage, cv::Mat& map, std::deque<int>& mapIndicesX,
	std::deque<int>& mapIndicesY, int low, int high)
{
	// 0 -> the pixel may be edge
	// 1 -> the pixel is not edge
	// 2 -> the pixel is edge
	map = cv::Mat::ones(NMSImage.size(), CV_8UC1);

	for (int i = 0; i < NMSImage.rows; ++i)
	{
		for (int j = 0; j < NMSImage.cols; ++j)
		{
			int m = NMSImage.ptr<uchar>(i)[j];				//nms -> CV_8UC1
			if (m > low)
			{
				if (m > high)
				{
					map.ptr<uchar>(i)[j] = 2;
					mapIndicesX.push_back(j);
					mapIndicesY.push_back(i);
				}
				else
					map.ptr<uchar>(i)[j] = 0;
			}
		}
	}
}


//双阈值滞后处理：根据队列中的像素坐标，进行8领域边缘点寻找，即在map中与2相连的0均认作为边缘点
void ImageDetectMethod::_hysteresis_thresholding(std::deque<int>& mapIndicesX, std::deque<int>& mapIndicesY, cv::Mat& map)
{
	while (!mapIndicesX.empty())
	{
		int r = mapIndicesY.back();
		int c = mapIndicesX.back();
		//获取到边缘点之后要将其弹出
		mapIndicesX.pop_back();
		mapIndicesY.pop_back();

		// top left
		if (map.ptr<uchar>(r - 1)[c - 1] == 0)
		{
			mapIndicesX.push_back(c - 1);
			mapIndicesY.push_back(r - 1);
			map.ptr<uchar>(r - 1)[c - 1] = 2;
		}
		// top
		if (map.ptr<uchar>(r - 1)[c] == 0)
		{
			mapIndicesX.push_back(c);
			mapIndicesY.push_back(r - 1);
			map.ptr<uchar>(r - 1)[c] = 2;
		}
		// top right
		if (map.ptr<uchar>(r - 1)[c + 1] == 0)
		{
			mapIndicesX.push_back(c + 1);
			mapIndicesY.push_back(r - 1);
			map.ptr<uchar>(r - 1)[c + 1] = 2;
		}
		// left
		if (map.ptr<uchar>(r)[c - 1] == 0)
		{
			mapIndicesX.push_back(c - 1);
			mapIndicesY.push_back(r);
			map.ptr<uchar>(r)[c - 1] = 2;
		}
		// right
		if (map.ptr<uchar>(r)[c + 1] == 0)
		{
			mapIndicesX.push_back(c + 1);
			mapIndicesY.push_back(r);
			map.ptr<uchar>(r)[c + 1] = 2;
		}
		// bottom left
		if (map.ptr<uchar>(r + 1)[c - 1] == 0)
		{
			mapIndicesX.push_back(c - 1);
			mapIndicesY.push_back(r + 1);
			map.ptr<uchar>(r + 1)[c - 1] = 2;
		}
		// bottom
		if (map.ptr<uchar>(r + 1)[c] == 0)
		{
			mapIndicesX.push_back(c);
			mapIndicesY.push_back(r + 1);
			map.ptr<uchar>(r + 1)[c] = 2;
		}
		// bottom right
		if (map.ptr<uchar>(r + 1)[c + 1] == 0)
		{
			mapIndicesX.push_back(c + 1);
			mapIndicesY.push_back(r + 1);
			map.ptr<uchar>(r + 1)[c + 1] = 2;
		}
	}
}
cv::Mat ImageDetectMethod::_get_canny_result(const cv::Mat& map)
{
	cv::Mat dst(map.rows - 2, map.cols - 2, CV_8UC1);
	for (int i = 0; i < dst.rows; i++) {
		for (int j = 0; j < dst.cols; j++) {
			dst.ptr<uchar>(i)[j] = (map.ptr<uchar>(i + 1)[j + 1] == 2 ? 255 : 0);
		}
	}
	return dst;
}
cv::Mat ImageDetectMethod::Adaptive_Canny(const cv::Mat& src, int apertureSize, bool L2gradient)
{
	CV_Assert(src.type() == CV_8UC1);
	CV_Assert(apertureSize == 3 || apertureSize == 5);

	cv::Mat dx, dy, magnitudes, angles;
	cv::Mat gaussianSrc = src.clone();
	_sobel_gradient(gaussianSrc, dx, dy, magnitudes, angles, apertureSize, L2gradient);

	//非极大值抑制计算高低阈值
	int low, high;
	cv::Mat NMSImage;
	_calculate_hysteresis_threshold_value(dx, dy, magnitudes, angles, NMSImage, low, high);

	cv::Mat map;
	std::deque<int> mapIndicesX, mapIndicesY;
	_non_maximum_suppression(NMSImage, map, mapIndicesX, mapIndicesY, low, high);

	_hysteresis_thresholding(mapIndicesX, mapIndicesY, map);
	cv::Mat dst = _get_canny_result(map);

	return dst;
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
bool ImageDetectMethod::TraceEdge(int start_y, int start_x, int y, int x,
	std::vector<cv::Point>* edge_vector, cv::Mat* cannyed_image_mat, cv::Mat* result_mat, bool* is_contour_closed)
{
	//对8邻域像素进行查询,搜索方案 直角边，斜对角边
// 	int xNum[8] = {1,1,0,-1,-1,-1,0,1};
// 	int yNum[8] = {0,1,1,1,0,-1,-1,-1};
	int xNum[8] = { 1, 0,-1, 0, 1,-1,-1,1 };
	int yNum[8] = { 0, 1, 0,-1, 1, 1,-1,-1 };
	// 	int xNum[8] = {1,-1,-1, 1, 1, 0,-1, 0};
	// 	int yNum[8] = {1, 1,-1,-1, 0, 1, 0,-1 };

	int yy, xx;

	for (int k = 0; k < 8; k++)
	{
		yy = y + yNum[k];
		xx = x + xNum[k];

		if (xx<1 || xx>cannyed_image_mat->cols - 1 || yy<1 || yy>cannyed_image_mat->rows - 1)
		{
			continue;
		}

		//
		if (*is_contour_closed == true)
		{
			return true;
		}
		//递归深度限制，不然会堆栈溢出
		if (edge_vector->size() > 2000)
		{
			return false;
		}

		uchar* _cannyed_image_mat_ptr = cannyed_image_mat->ptr<uchar>(yy);
		uchar* _result_mat_ptr = result_mat->ptr<uchar>(yy);

		if (_cannyed_image_mat_ptr[xx] > 0 && _result_mat_ptr[xx] == 0)
		{
			//该点设为边界点
			_result_mat_ptr[xx] = 255;
			edge_vector->push_back(cv::Point(xx, yy));

			//判断是否不是中间有跳跃
			if (edge_vector->size() > 3)
			{
				if (abs(edge_vector->at(edge_vector->size() - 1).x - edge_vector->at(edge_vector->size() - 2).x) +
					abs(edge_vector->at(edge_vector->size() - 1).y - edge_vector->at(edge_vector->size() - 2).y) > 2)
				{
					return false;
				}
			}

			//以该点为中心再进行跟踪
			if (TraceEdge(start_y, start_x, yy, xx, edge_vector, cannyed_image_mat, result_mat, is_contour_closed) == false)
			{
				return false;
			}
		}
	}

	//判断轮廓闭合
	if (edge_vector->size() > 2 && sqrt(double(edge_vector->at(edge_vector->size() - 1).x - start_x) * double(edge_vector->at(edge_vector->size() - 1).x - start_x)
		+ double(edge_vector->at(edge_vector->size() - 1).y - start_y) * double(edge_vector->at(edge_vector->size() - 1).y - start_y)) < 1.5)
	{
		*is_contour_closed = true;
	}

	return true;
}

bool ImageDetectMethod::SelfCannyMethod(const cv::Mat& ori_image_mat, cv::Mat& output_image_mat,
	float canny_low_thresh, float canny_high_thresh, std::vector<std::vector<cv::Point>>& edge_point_list)
{
	cv::Mat src;
	ori_image_mat.convertTo(src, CV_32F);

	int mat_rows = ori_image_mat.rows;
	int mat_cols = ori_image_mat.cols;

	cv::Mat M_gradient = cv::Mat::zeros(mat_rows, mat_cols, CV_32F);
	cv::Mat sita_mat = cv::Mat::zeros(mat_rows, mat_cols, CV_32F);

	if (canny_low_thresh > canny_high_thresh)
		std::swap(canny_low_thresh, canny_high_thresh);
	int low = cvFloor(canny_low_thresh);
	int high = cvFloor(canny_high_thresh);

	cv::Mat dx(src.rows, src.cols, CV_32F);
	cv::Mat dy(src.rows, src.cols, CV_32F);

	cv::Sobel(src, dx, CV_32F, 1, 0, 3, 1, 0, cv::BORDER_REPLICATE);
	cv::Sobel(src, dy, CV_32F, 0, 1, 3, 1, 0, cv::BORDER_REPLICATE);

	for (int i = 1; i < src.rows - 1; i++)
	{
		float* _norm = M_gradient.ptr<float>(i);
		float* _sita = sita_mat.ptr<float>(i);
		float* _dx = dx.ptr<float>(i);
		float* _dy = dy.ptr<float>(i);

		for (int j = 1; j < src.cols - 1; j++)
		{
			//_norm[j] = std::abs(_dx[j]) + std::abs(_dy[j]);
			_norm[j] = sqrt(_dx[j] * _dx[j] + _dy[j] * _dy[j]);
			_sita[j] = atan2(_dy[j], _dx[j]);
		}
	}

	//Reconstruct3DMethod::WirteMatToFile(M_gradient);

	cv::Mat NMS_mat = cv::Mat::zeros(mat_rows, mat_cols, CV_32F);

	for (int i = 1; i < mat_rows - 1; i++)
	{
		float* _NMS = NMS_mat.ptr<float>(i);
		float* _sita = sita_mat.ptr<float>(i);
		float* _norm = M_gradient.ptr<float>(i);
		float* _norm_up = M_gradient.ptr<float>(i - 1);
		float* _norm_down = M_gradient.ptr<float>(i + 1);
		for (int j = 1; j < mat_cols - 1; j++)
		{
			float m1, m2;                //梯度方向相邻的两个梯度值

			if (abs(_sita[j]) <= M_PI / 8 || abs(_sita[j]) >= 7 * M_PI / 8)
			{
				m1 = _norm[j + 1];
				m2 = _norm[j - 1];
			}
			else if ((_sita[j] >= M_PI / 8 && _sita[j] < 3 * M_PI / 8) || (_sita[j] > -7 * M_PI / 8 && _sita[j] < -5 * M_PI / 8))
			{
				m1 = _norm_down[j + 1];
				m2 = _norm_up[j - 1];
			}
			else if (abs(_sita[j]) >= 3 * M_PI / 8 && abs(_sita[j]) <= 5 * M_PI / 8)
			{
				m1 = _norm_down[j];
				m2 = _norm_up[j];
			}
			else
			{
				m1 = _norm_up[j + 1];
				m2 = _norm_down[j - 1];
			}

			if (_norm[j] >= m1 && _norm[j] >= m2)
			{
				_NMS[j] = float(_norm[j]);
			}
		}
	}

	/////Reconstruct3DMethod::WirteMatToFile(NMS_mat);

	///////////双阈值,指定阈值
	// 	Mat low_threshold_mat;
	// 	Mat high_threhold_mat;
	// 	threshold(NMS_mat, low_threshold_mat, low, 255, THRESH_BINARY);
	// 	threshold(NMS_mat, high_threhold_mat, high, 255, THRESH_BINARY);

	/////canny自适应阈值
	////otsu法，先将非极大值抑制后梯度的值进行，归一化（整数化），分成64等级，到高阈值th2，th1 = 0.5th2；
	int N = 64;
	cv::Mat norm_NMS_mat = cv::Mat(mat_rows, mat_cols, CV_8U);
	double max_value;
	minMaxLoc(NMS_mat, NULL, &max_value);
	double per_d = (max_value + 0.000001) / N;

	for (int i = 0; i < mat_rows; i++)
	{
		float* _NMS_mat = NMS_mat.ptr<float>(i);
		uchar* _norm_NMS_mat = norm_NMS_mat.ptr<uchar>(i);
		for (int j = 0; j < mat_cols; j++)
		{
			_norm_NMS_mat[j] = floor(_NMS_mat[j] / per_d);
		}
	}

	double max_level = OTSUForCannyAdaptiveHighTh(norm_NMS_mat, N);
	double high_th = max_level * per_d;
	double low_th = 0.5 * high_th;
	cv::Mat low_threshold_mat;
	cv::Mat high_threhold_mat;
	threshold(NMS_mat, low_threshold_mat, low_th, 255, cv::THRESH_BINARY);
	threshold(NMS_mat, high_threhold_mat, high_th, 255, cv::THRESH_BINARY);

	///////////
	// 	namedWindow("canny",CV_WINDOW_FREERATIO);
	// 	imshow("canny",low_threshold_mat);
	// 	imwrite("E:\\cannylow.bmp",low_threshold_mat);
	// 	waitKey(0);
	//
	// 	namedWindow("canny",CV_WINDOW_FREERATIO);
	// 	imshow("canny",high_threhold_mat);
	// 	imwrite("E:\\cannyhigh.bmp",high_threhold_mat);
	// 	waitKey(0);
	//
	// 	int ele[3][3] ={{1,1,1},{1,1,1},{1,1,1}};
	// 	Mat element = Mat::ones(5,5,CV_8U);
	// 	Mat low_threshold_mat_changed;
	//
	// 	morphologyEx(low_threshold_mat,low_threshold_mat,MORPH_CLOSE,element );
	// 	///膨胀
	// 	//dilate(low_threshold_mat,low_threshold_mat_changed,element);
	// 	///腐蚀
	// 	//erode（low_threshold_mat,low_threshold_mat_changed,element);
	//
	// 	namedWindow("canny",CV_WINDOW_FREERATIO);
	// 	imshow("canny",low_threshold_mat);
	// 	imwrite("E:\\cannylowd.bmp",low_threshold_mat);
	// 	waitKey(0);

	///////////////////////

	cv::Mat search_sign = cv::Mat::zeros(mat_rows, mat_cols, CV_8U);

	for (int i = 0; i < mat_rows; i++)
	{
		for (int j = 0; j < mat_cols; j++)
		{
			if (search_sign.at<uchar>(i, j) == 0 && high_threhold_mat.at<float>(i, j) > 0)
			{
				std::vector<cv::Point> edge_point;
				bool is_contour_closed = false;
				search_sign.at<uchar>(i, j) = 255;
				edge_point.push_back(cv::Point(j, i));

				FindCannyEdge(i, j, i, j, &low_threshold_mat, &high_threhold_mat, &search_sign, &edge_point, &is_contour_closed);

				if (is_contour_closed == true)
				{
					edge_point_list.push_back(edge_point);
				}
			}
		}
	}

	output_image_mat = search_sign;
	// 	gray_gradient_value = M_gradient;
	// 	gray_gradient_sita = sita_mat;

	return true;
}

bool ImageDetectMethod::FindCannyEdge(int start_y, int start_x, int y, int x,
	cv::Mat* low_threshold_mat, cv::Mat* high_threshold_mat, cv::Mat* search_sign_mat,
	std::vector<cv::Point>* edge_point_vector, bool* is_contour_closed)
{
	//对8邻域像素进行查询,搜索方案 直角边，斜对角边
	// 	int xNum[8] = {1,1,0,-1,-1,-1,0,1};
	// 	int yNum[8] = {0,1,1,1,0,-1,-1,-1};
	int xNum[8] = { 1, 0,-1, 0, 1,-1,-1,1 };
	int yNum[8] = { 0, 1, 0,-1, 1, 1,-1,-1 };
	// 	int xNum[8] = {1,-1,-1, 1, 1, 0,-1, 0};
	// 	int yNum[8] = {1, 1,-1,-1, 0, 1, 0,-1 };

	int i, j;

	bool is_found_next_point = false;
	for (int k = 0; k < 8; k++)
	{
		i = y + yNum[k];
		j = x + xNum[k];

		if (j<0 || j>search_sign_mat->cols - 1 || i<0 || i>search_sign_mat->rows - 1)
		{
			continue;
		}
		//判断轮廓闭合
		if (*is_contour_closed == true)
		{
			return true;
		}
		if (edge_point_vector->size() > 2000)
		{
			return false;
		}

		float* _high_threshold_mat_ptr = high_threshold_mat->ptr<float>(i);
		uchar* _search_sign_mat_ptr = search_sign_mat->ptr<uchar>(i);

		if (_high_threshold_mat_ptr[j] > 0 && _search_sign_mat_ptr[j] == 0)
		{
			if (edge_point_vector->size() > 3)
			{
				if (abs(edge_point_vector->at(edge_point_vector->size() - 1).x - edge_point_vector->at(edge_point_vector->size() - 2).x) +
					abs(edge_point_vector->at(edge_point_vector->size() - 1).y - edge_point_vector->at(edge_point_vector->size() - 2).y) > 2)
				{
					return false;
				}
			}

			is_found_next_point = true;
			_search_sign_mat_ptr[j] = 255;
			edge_point_vector->push_back(cv::Point(j, i));

			if (FindCannyEdge(start_y, start_x, i, j, low_threshold_mat, high_threshold_mat, search_sign_mat, edge_point_vector, is_contour_closed) == false)
			{
				return false;
			}
		}
	}

	if (is_found_next_point == false)
	{
		for (int k = 0; k < 8; k++)
		{
			i = y + yNum[k];
			j = x + xNum[k];

			if (j<0 || j>search_sign_mat->cols - 1 || i<0 || i>search_sign_mat->rows - 1)
			{
				continue;
			}
			//判断轮廓闭合
			if (*is_contour_closed == true)
			{
				return true;
			}
			if (edge_point_vector->size() > 2000)
			{
				return false;
			}

			float* _low_threshold_mat_ptr = low_threshold_mat->ptr<float>(i);
			uchar* _search_sign_mat_ptr = search_sign_mat->ptr<uchar>(i);

			if (_search_sign_mat_ptr[j] > 0)
			{
				continue;
			}

			if (_low_threshold_mat_ptr[j] > 0 && _search_sign_mat_ptr[j] == 0)
			{
				if (edge_point_vector->size() > 3)
				{
					if (abs(edge_point_vector->at(edge_point_vector->size() - 1).x - edge_point_vector->at(edge_point_vector->size() - 2).x) +
						abs(edge_point_vector->at(edge_point_vector->size() - 1).y - edge_point_vector->at(edge_point_vector->size() - 2).y) > 2)
					{
						return false;
					}
				}

				if (edge_point_vector->size() > 2000)
				{
					break;
				}

				is_found_next_point = true;
				_search_sign_mat_ptr[j] = 255;
				edge_point_vector->push_back(cv::Point(j, i));

				if (FindCannyEdge(start_y, start_x, i, j, low_threshold_mat, high_threshold_mat, search_sign_mat, edge_point_vector, is_contour_closed) == false)
				{
					return false;
				}
			}
		}
	}

	//判断轮廓闭合
	if (edge_point_vector->size() > 2 && sqrt(double(edge_point_vector->at(edge_point_vector->size() - 1).x - start_x) * double(edge_point_vector->at(edge_point_vector->size() - 1).x - start_x)
		+ double(edge_point_vector->at(edge_point_vector->size() - 1).y - start_y) * double(edge_point_vector->at(edge_point_vector->size() - 1).y - start_y)) < 1.5)
	{
		*is_contour_closed = true;
	}

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
bool ImageDetectMethod::EllipseGrayJudgeForPointCSI2(const cv::Mat& image_mat, float center_x, float center_y,
	float ellipse_a, float ellipse_b, float angle_in_pi, float ratio_k)
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

	QList<float> all_value_list_0;
	QList<float> all_value_list_01;
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
				all_value_list.append(float(_image_mat_ptr[j]));
				all_value_list_0.append(-1);
				all_value_list_01.append(-1);
				//std::cout << 0 << "\t" << i << "\t" << j << "\t" << float(_image_mat_ptr[j]) << "\n";
			}
			else if (tr_x * tr_x / int(ellipse_a) / int(ellipse_a) + tr_y * tr_y / int(ellipse_b) / int(ellipse_b) < 1)
			{
				all_value_list.append(float(_image_mat_ptr[j]));
				all_value_list_0.append(-1);
				all_value_list_01.append(1);
			}
			else
			{
				all_value_list.append(float(_image_mat_ptr[j]));
				all_value_list_0.append(1);
				all_value_list_01.append(1);
			}
		}
	}
	double mean_value_ori = MeanValue(all_value_list);
	double mean_value_00 = MeanValue(all_value_list_0);
	double mean_value_01 = MeanValue(all_value_list_01);
	double std_value_ori = 0;
	double std_value_00 = 0;
	double std_value_01 = 0;
	for (int ii = 0; ii < all_value_list.size(); ii++)
	{
		std_value_ori += pow(all_value_list[ii] - mean_value_ori, 2);
		std_value_00 += pow(all_value_list_0[ii] - mean_value_00, 2);
		std_value_01 += pow(all_value_list_01[ii] - mean_value_01, 2);
	}
	std_value_ori = sqrt(std_value_ori / (double)all_value_list.size());
	std_value_00 = sqrt(std_value_00 / (double)all_value_list_0.size());
	std_value_01 = sqrt(std_value_01 / (double)all_value_list_01.size());
	double zncc_00 = 0;
	double zncc_01 = 0;

	for (int ii = 0; ii < all_value_list.size(); ii++)
	{
		zncc_00 += (all_value_list[ii] - mean_value_ori) * (all_value_list_0[ii] - mean_value_00) / std_value_ori / std_value_00;
		zncc_01 += (all_value_list[ii] - mean_value_ori) * (all_value_list_01[ii] - mean_value_01) / std_value_ori / std_value_01;
	}
	zncc_00 /= (double)all_value_list.size();
	zncc_01 /= (double)all_value_list.size();
	if (abs(zncc_01) > abs(zncc_01))
	{
		return true;
	}
	return false;
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
		//int mid_value = AverageOfList(gray_list);
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

bool ImageDetectMethod::FindSubPixelPosOfCircleCenter20140210(const cv::Mat& image_mat, float center_x, float center_y, float ellipse_a, float ellipse_b,
	float angle_in_pi, const std::vector<cv::Point>& contour_points,
	float& sub_pixel_center_x, float& sub_pixel_center_y,
	std::vector<cv::Point2f>* subpixel_edge_points /*= NULL*/,
	SubPixelPosMethod subPixel_method /*= NoSubPixel_Match*/,
	MarkPointColorType color_type/* = BlackDownWhiteUp*/)
{
	if (subPixel_method == NoSubPixel_Match)
	{
		sub_pixel_center_x = center_x;
		sub_pixel_center_y = center_y;

		if (subpixel_edge_points != NULL)
		{
			for (int i = 0; i < contour_points.size(); i++)
			{
				subpixel_edge_points->push_back(cv::Point2f(contour_points.at(i).x, contour_points.at(i).y));
			}
		}
		return true;
	}
	else
		if (subPixel_method == Binary_Centroid || subPixel_method == Gray_Centroid || subPixel_method == Squared_Gray_Centroid)
		{
			//重心法,局部灰度阈值需要确认
			QRect sub_rect = GetEllipseROIRect(image_mat, center_x, center_y, ellipse_a, ellipse_b, angle_in_pi);
			cv::Mat sub_mat = image_mat.operator ()(cv::Rect(sub_rect.x(), sub_rect.y(), sub_rect.width() + 1, sub_rect.height() + 1));

			if (color_type == Uncertainty)
			{
				color_type = JudgeTargetColorType(sub_mat, center_x - sub_rect.x(), center_y - sub_rect.y(),
					ellipse_a, ellipse_b, angle_in_pi);
			}

			//阈值选择方法，还有固定阈值 ，不加阈值0,255,暂时不知道哪种阈值方法较好，需要测试
			float gray_threshold;
			// 		switch(color_type)                //固定阈值值
			// 		{
			// 		case BlackDownWhiteUp:
			// 			gray_threshold = 0;
			// 			break;
			// 		case WhiteDownBlackUp:
			// 			gray_threshold =255;
			// 			break;
			// 		}
			gray_threshold = DICOTSU20140215(sub_mat);   //otsu算法确定阈值
			/*gray_threshold = CalThresholdInSubsetMat(sub_mat,center_x - sub_rect.x(),center_y -sub_rect.y(),
			   ellipse_a,ellipse_b,angle_in_pi);
		   gray_threshold= CalThresholdInSubsetMat2(image_mat,contour_points);*/

#pragma region 计算圆心,选择二值化灰度重心，灰度重心法，灰度平方重心法
		   //计算圆心,选择二值化灰度重心，灰度重心法，灰度平方重心法
   // 		float sub_x,sub_y;
   // 		if (CalCentriodBySubsetMat(sub_mat,gray_threshold,sub_x,sub_y,color_type,subPixel_method))
   // 		{
   // 			sub_pixel_center_x = sub_x + sub_rect.x();
   // 			sub_pixel_center_y = sub_y + sub_rect.y();
   //
   //
   // 			if (subpixel_edge_points !=NULL)
   // 			{
   // 				for (int i=0;i<contour_points.size();i++)
   // 				{
   // 					subpixel_edge_points->push_back(Point2f(contour_points.at(i).x,contour_points.at(i).y));
   // 				}
   // 			}
   //
   // 			return true;
   // 		}else
   // 			return false;
#pragma endregion
		//**************多阈值灰度形心法、灰度重心法，灰度平方重心法
			int k = 0;
			int d_gray = 6;
			float sub_x, sub_y;
			sub_x = sub_y = 0;
			for (int i = -k; i <= k; i++)
			{
				float threshold = gray_threshold + i * d_gray;
				float part_sub_x, part_sub_y;
				CalCentriodBySubsetMat(sub_mat, threshold, part_sub_x, part_sub_y, color_type, subPixel_method);
				sub_x += part_sub_x;
				sub_y += part_sub_y;
			}
			sub_x /= 2 * k + 1;
			sub_y /= 2 * k + 1;

			sub_pixel_center_x = sub_x + sub_rect.x();
			sub_pixel_center_y = sub_y + sub_rect.y();

			if (subpixel_edge_points != NULL)
			{
				for (int i = 0; i < contour_points.size(); i++)
				{
					subpixel_edge_points->push_back(cv::Point2f(contour_points.at(i).x, contour_points.at(i).y));
				}
			}

			return true;
		}
		else
			if (subPixel_method == Interpolation_Ellipse_Match)
			{
				//<<亚像素边缘定位算法的稳定性分析>> 田原嫄,2010
				// x和y方向插值，效果不好，x-1,x,x+1的梯度fu值x处最大时较好
				//二次多项式插值，再椭圆拟合
				std::vector<cv::Point2f> new_contour_points;

				for (int n = 0; n < contour_points.size(); n++)
				{
					int y = contour_points[n].y;
					int x = contour_points[n].x;

					//    0  1  2
					//    3  4  5
					//    6  7  8
					float gray_gradient[9];
					float gradient_sita;
					CalGrayGradientBySobel(image_mat, y - 1, x, &gray_gradient[1]);
					CalGrayGradientBySobel(image_mat, y, x - 1, &gray_gradient[3]);
					CalGrayGradientBySobel(image_mat, y, x, &gray_gradient[4], &gradient_sita);
					CalGrayGradientBySobel(image_mat, y, x + 1, &gray_gradient[5]);
					CalGrayGradientBySobel(image_mat, y + 1, x, &gray_gradient[7]);

					float dx, dy;
					float dx_den = (gray_gradient[3] - 2 * gray_gradient[4] + gray_gradient[5]);
					float dy_den = (gray_gradient[1] - 2 * gray_gradient[4] + gray_gradient[7]);

					if (dx_den == 0 || dy_den == 0)
					{
						continue;
					}

					if (gray_gradient[4] > gray_gradient[1] && gray_gradient[4] > gray_gradient[7]
						&& gray_gradient[4] > gray_gradient[3] && gray_gradient[4] > gray_gradient[5])
					{
						dx = (gray_gradient[3] - gray_gradient[5]) / 2 / dx_den;
						dy = (gray_gradient[1] - gray_gradient[7]) / 2 / dy_den;
					}
					else
						continue;

					new_contour_points.push_back(cv::Point2f(x + dx, y + dy));
				}

				cv::RotatedRect new_rRect = fitEllipse(new_contour_points);    //fitEllipse只接受float和int类型

				//根据误差剔除部分点
				ReduceBadEllipseFitPoints(new_contour_points, new_rRect.center.x, new_rRect.center.y,
					new_rRect.size.width * 0.5, new_rRect.size.height * 0.5, new_rRect.angle * M_PI / 180);
				if (new_contour_points.size() < 6)
				{
					return false;
					// 			sub_pixel_center_x = 0;
					// 			sub_pixel_center_y = 0;
					// 			return true;
				}
				new_rRect = fitEllipse(new_contour_points);

				sub_pixel_center_x = new_rRect.center.x;
				sub_pixel_center_y = new_rRect.center.y;

				if (subpixel_edge_points != NULL)
				{
					*subpixel_edge_points = new_contour_points;
				}

				return true;
			}
			else
				if (subPixel_method == Interpolation_Rotaion_Ellipse_Match)
				{
					//二次多项式插值，坐标转换后，只插值梯度方向，后进行坐标转换得到dx，dy
					std::vector<cv::Point2f> new_contour_points;

					for (int n = 0; n < contour_points.size(); n++)
					{
						int y = contour_points[n].y;
						int x = contour_points[n].x;

						//    0  1  2
						//    3  4  5
						//    6  7  8
						float gray_gradient[9];
						float gradient_sita;
						CalGrayGradientBySobel(image_mat, y - 1, x, &gray_gradient[1]);
						CalGrayGradientBySobel(image_mat, y, x - 1, &gray_gradient[3]);
						CalGrayGradientBySobel(image_mat, y, x, &gray_gradient[4], &gradient_sita);
						CalGrayGradientBySobel(image_mat, y, x + 1, &gray_gradient[5]);
						CalGrayGradientBySobel(image_mat, y + 1, x, &gray_gradient[7]);

						//    2          | -->x
						// 1  0  3       |
						//    4           y
						float boder_gradient_value[5];
						boder_gradient_value[0] = gray_gradient[4];
						//插值求算

						//简单求算，与canny算子非极大值抑制
						if ((gradient_sita > -M_PI / 8 && gradient_sita < -7 * M_PI / 8) || (gradient_sita >= 7 * M_PI / 8 && gradient_sita <= M_PI))
						{
							boder_gradient_value[1] = gray_gradient[5];
							boder_gradient_value[2] = gray_gradient[3];
							boder_gradient_value[3] = gray_gradient[7];
							boder_gradient_value[4] = gray_gradient[1];
						}
						else if (gradient_sita >= -7 * M_PI / 8 && gradient_sita < -5 * M_PI / 8)
						{
							boder_gradient_value[1] = gray_gradient[8];
							boder_gradient_value[2] = gray_gradient[0];
							boder_gradient_value[3] = gray_gradient[6];
							boder_gradient_value[4] = gray_gradient[2];
						}
						else if (gradient_sita >= -5 * M_PI / 8 && gradient_sita < -3 * M_PI / 8)
						{
							boder_gradient_value[1] = gray_gradient[7];
							boder_gradient_value[2] = gray_gradient[1];
							boder_gradient_value[3] = gray_gradient[3];
							boder_gradient_value[4] = gray_gradient[5];
						}
						else if (gradient_sita >= -3 * M_PI / 8 && gradient_sita < -M_PI / 8)
						{
							boder_gradient_value[1] = gray_gradient[6];
							boder_gradient_value[2] = gray_gradient[2];
							boder_gradient_value[3] = gray_gradient[0];
							boder_gradient_value[4] = gray_gradient[8];
						}
						else if (gradient_sita >= -M_PI / 8 && gradient_sita < M_PI / 8)
						{
							boder_gradient_value[1] = gray_gradient[3];
							boder_gradient_value[2] = gray_gradient[5];
							boder_gradient_value[3] = gray_gradient[1];
							boder_gradient_value[4] = gray_gradient[7];
						}
						else if (gradient_sita >= M_PI / 8 && gradient_sita < 3 * M_PI / 8)
						{
							boder_gradient_value[1] = gray_gradient[0];
							boder_gradient_value[2] = gray_gradient[8];
							boder_gradient_value[3] = gray_gradient[2];
							boder_gradient_value[4] = gray_gradient[6];
						}
						else if (gradient_sita >= 3 * M_PI / 8 && gradient_sita < 5 * M_PI / 8)
						{
							boder_gradient_value[1] = gray_gradient[1];
							boder_gradient_value[2] = gray_gradient[7];
							boder_gradient_value[3] = gray_gradient[5];
							boder_gradient_value[4] = gray_gradient[3];
						}
						else
						{
							boder_gradient_value[1] = gray_gradient[2];
							boder_gradient_value[2] = gray_gradient[6];
							boder_gradient_value[3] = gray_gradient[8];
							boder_gradient_value[4] = gray_gradient[0];
						}

						//将xy坐标系顺时针旋转sita角，将x正向指向梯度正向（及由圆心向外发射方向）
						float dx, dy;
						float dx_den = boder_gradient_value[1] - 2 * boder_gradient_value[0] + boder_gradient_value[3];
						float dy_den = boder_gradient_value[2] - 2 * boder_gradient_value[0] + boder_gradient_value[4];

						if (dx_den == 0 || dy_den == 0)
						{
							continue;
						}

						if (gray_gradient[4] > gray_gradient[1] && gray_gradient[4] > gray_gradient[7]
							&& gray_gradient[4] > gray_gradient[3] && gray_gradient[4] > gray_gradient[8])
						{
							dx = (boder_gradient_value[1] - boder_gradient_value[3]) / 2 / dx_den;
							dy = (boder_gradient_value[2] - boder_gradient_value[4]) / 2 / dy_den;
						}
						else
							continue;

						//坐标转换为正常坐标
						float dx_re = dx * cos(gradient_sita) + dy * sin(gradient_sita);
						float dy_re = -dx * sin(gradient_sita) + dy * cos(gradient_sita);

						new_contour_points.push_back(cv::Point2f(x + dx_re, y + dy_re));
					}
					cv::RotatedRect new_rRect = fitEllipse(new_contour_points);    //fitEllipse只接受float和int类型

					//根据误差剔除部分点
					ReduceBadEllipseFitPoints(new_contour_points, new_rRect.center.x, new_rRect.center.y,
						new_rRect.size.width * 0.5, new_rRect.size.height * 0.5, new_rRect.angle * M_PI / 180);
					if (new_contour_points.size() < 6)
					{
						return false;
					}
					new_rRect = fitEllipse(new_contour_points);

					sub_pixel_center_x = new_rRect.center.x;
					sub_pixel_center_y = new_rRect.center.y;

					return true;
				}
				else
					if (subPixel_method == Gauss_surface_fit_Ellipse_Match)
					{
						//高斯曲面拟合法
						std::vector<cv::Point2f> new_contour_points;

						for (int n = 0; n < contour_points.size(); n++)
						{
							int y = contour_points[n].y;
							int x = contour_points[n].x;

							cv::Mat A = cv::Mat(9, 5, CV_64F);
							cv::Mat B = cv::Mat(9, 1, CV_64F);
							for (int i = -1; i < 2; i++)
							{
								for (int j = -1; j < 2; j++)
								{
									A.at<double>((i + 1) * 3 + j + 1, 0) = j * j * double(image_mat.at<uchar>(y + i, x + j));
									A.at<double>((i + 1) * 3 + j + 1, 1) = i * i * double(image_mat.at<uchar>(y + i, x + j));
									A.at<double>((i + 1) * 3 + j + 1, 2) = j * double(image_mat.at<uchar>(y + i, x + j));
									A.at<double>((i + 1) * 3 + j + 1, 3) = i * double(image_mat.at<uchar>(y + i, x + j));
									A.at<double>((i + 1) * 3 + j + 1, 4) = double(image_mat.at<uchar>(y + i, x + j));

									B.at<double>((i + 1) * 3 + j + 1, 0) = double(image_mat.at<uchar>(y + i, x + j)) * log(double(image_mat.at<uchar>(y + i, x + j)));
								}
							}

							cv::Mat X;
							solve(A.t() * A, A.t() * B, X);

							double dx = -X.at<double>(2, 0) / 2 / X.at<double>(0, 0);
							double dy = -X.at<double>(3, 0) / 2 / X.at<double>(1, 0);

							double a[9][5];
							for (int i = 0; i < 9; i++)
							{
								for (int j = 0; j < 5; j++)
								{
									a[i][j] = A.at<double>(i, j);
								}
							}
							double b[9];
							for (int i = 0; i < 9; i++)
							{
								b[i] = B.at<double>(i, 0);
							}
							double xx[5];
							for (int i = 0; i < 5; i++)
							{
								xx[i] = X.at<double>(i, 0);
							}

							if (_isnan(dx) || _isnan(dy))
							{
								continue;
							}

							new_contour_points.push_back(cv::Point2f(x + dx, y + dy));
						}
						cv::RotatedRect new_rRect = fitEllipse(new_contour_points);    //fitEllipse只接受float和int类型

						sub_pixel_center_x = new_rRect.center.x;
						sub_pixel_center_y = new_rRect.center.y;

						return true;
					}
					else
						if (subPixel_method == Surface_fit_Ellipse_Match)
						{
							// <<圆形标志点的亚像素定位及其应用>>  殷永凯，2008
							//曲面拟合法，对灰度进行曲面拟合，拟合公式 f(x,y)= k1+k2x+k3y+k4x^2+k5xy+k6y^2+k7x^3+k8x^2y+k9xy^2+k10y^3
							//沿梯度方向求取二阶导为零值
							std::vector<cv::Point2f> new_contour_points;

							for (int n = 0; n < contour_points.size(); n++)
							{
								int y = contour_points[n].y;
								int x = contour_points[n].x;

								if (x - 3 < 0 || x + 3 > image_mat.cols - 1 || y - 3 < 0 || y + 3 > image_mat.rows - 1)
								{
									continue;
								}

								int half_size = 2;
								int total_pixel = (2 * half_size + 1) * (2 * half_size + 1);
								cv::Mat A = cv::Mat(total_pixel, 10, CV_32F);
								cv::Mat B = cv::Mat(total_pixel, 1, CV_32F);

								for (int i = -half_size; i < half_size + 1; i++)
								{
									const uchar* _image_mat = image_mat.ptr<uchar>(i + y);
									for (int j = -half_size; j < half_size + 1; j++)
									{
										float* _A = A.ptr<float>((i + 2) * 5 + j + 2);
										float* _B = B.ptr<float>((i + 2) * 5 + j + 2);

										_A[0] = 1;
										_A[1] = j;
										_A[2] = i;
										_A[3] = j * j;
										_A[4] = j * i;
										_A[5] = i * i;
										_A[6] = j * j * j;
										_A[7] = j * j * i;
										_A[8] = j * i * i;
										_A[9] = i * i * i;

										_B[0] = float(_image_mat[j + x]);
									}
								}

								//X 10*1: [k1,k2,k3,k4,k5,k6,k7,k8,k9,k10]'
								cv::Mat X;
								solve(A.t() * A, A.t() * B, X);

								float gray_gradient, gradient_sita;
								CalGrayGradientBySobel(image_mat, y, x, &gray_gradient, &gradient_sita);

								cv::Mat XT = X.t();
								float* _XT_ptr = XT.ptr<float>(0);
								//f(x,y)二阶导=0，得到a*r+b=0,求r, d = sita
								//a = 6*(k7*sind^3 +k8*sind^2*cosd + k9sindcosd^2 +k10*cosd^3)
								// b= 2*(k4*sind^2 + k5*sind*cosd + k6*cosd^2)
								float a = 6 * (_XT_ptr[6] * sin(gradient_sita) * sin(gradient_sita) * sin(gradient_sita)
									+ _XT_ptr[7] * sin(gradient_sita) * sin(gradient_sita) * cos(gradient_sita)
									+ _XT_ptr[8] * sin(gradient_sita) * sin(gradient_sita) * cos(gradient_sita)
									+ _XT_ptr[9] * sin(gradient_sita) * sin(gradient_sita) * cos(gradient_sita));
								float b = 2 * (_XT_ptr[3] * sin(gradient_sita) * sin(gradient_sita)
									+ _XT_ptr[4] * sin(gradient_sita) * cos(gradient_sita)
									+ _XT_ptr[5] * cos(gradient_sita) * cos(gradient_sita));

								float r = -b / a;

								new_contour_points.push_back(cv::Point2f(x + r * cos(gradient_sita), y + r * sin(gradient_sita)));
							}

							cv::RotatedRect new_rRect = fitEllipse(new_contour_points);    //fitEllipse只接受float和int类型

							sub_pixel_center_x = new_rRect.center.x;
							sub_pixel_center_y = new_rRect.center.y;

							return true;
						}
						else
							if (subPixel_method == Gauss_Curve_Fit)
							{
								//<<光学测量中椭圆圆心定位算法研究>>  张虎， 2008
								/////沿梯度方向对灰度梯度幅值进行高斯曲线拟合
								//////先插值出梯度方向的灰度值，再用
								std::vector<cv::Point2f> new_contour_points;

								for (int n = 0; n < contour_points.size(); n++)
								{
									int y = contour_points[n].y;
									int x = contour_points[n].x;

									if (x - 3 < 0 || x + 3 > image_mat.cols - 1 || y - 3 < 0 || y + 3 > image_mat.rows - 1)
									{
										continue;
									}

									if (x - center_x == 0)
									{
										continue;
									}

									float sita = atan2(y - center_y, x - center_x);
									float k = tan(sita);
									//float k = (y-center_y)/(x-center_x);   //梯度方向
									// 			double sita = gray_gradient_sita.at<double>(y,x);
									// 			double k = tan(sita);
									//double f_a,f_b,f_c,f_d;
									float abs_k = abs(k);

									float positive_gray[3];
									float negative_gray[3];

									int invt = 1;
									if (k > 0)
									{
										invt = -1;
									}
									if (abs_k < 1)
									{
										for (int i = 1; i < 4; i++)
										{
											int a = int(i * abs_k) / 1;
											float lamd2 = i * abs_k - a;
											float lamd1 = 1 - lamd2;

											positive_gray[i - 1] = lamd1 * image_mat.at<uchar>(y - a, x + invt * i) + lamd2 * image_mat.at<uchar>(y - a - 1, x + invt * i);
											negative_gray[i - 1] = lamd1 * image_mat.at<uchar>(y + a, x - invt * i) + lamd2 * image_mat.at<uchar>(y + a + 1, x - invt * i);
										}
									}
									else
									{
										for (int i = 1; i < 4; i++)
										{
											int a = int(i / abs_k) / 1;
											float lamd2 = i / abs_k - a;
											float lamd1 = 1 - lamd2;

											positive_gray[i - 1] = lamd1 * image_mat.at<uchar>(y - i, x + invt * a) + lamd2 * image_mat.at<uchar>(y - i, x + invt * (a + 1));
											negative_gray[i - 1] = lamd1 * image_mat.at<uchar>(y + i, x - invt * a) + lamd2 * image_mat.at<uchar>(y + i, x - invt * (a + 1));
										}
									}

									float gray_value[7];    //沿灰度梯度方向 ，向外
									if (y < center_y)
									{
										gray_value[0] = negative_gray[2];
										gray_value[1] = negative_gray[1];
										gray_value[2] = negative_gray[0];
										gray_value[3] = image_mat.at<uchar>(y, x);
										gray_value[4] = positive_gray[0];
										gray_value[5] = positive_gray[1];
										gray_value[6] = positive_gray[2];
									}
									else
									{
										gray_value[0] = positive_gray[2];
										gray_value[1] = positive_gray[1];
										gray_value[2] = positive_gray[0];
										gray_value[3] = image_mat.at<uchar>(y, x);
										gray_value[4] = negative_gray[0];
										gray_value[5] = negative_gray[1];
										gray_value[6] = negative_gray[2];
									}

									//计算差分
									float f_difference[5];
									for (int i = 0; i < 5; i++)
									{
										f_difference[i] = abs(gray_value[i] - gray_value[i + 1]) / 2 + abs(gray_value[i + 1] - gray_value[i + 2]) / 2; //向前插值和向后插值
										//f_difference[i] = abs(gray_value[i]-gray_value[i+2]);
									}

									//高斯曲线拟合
									float delta = (0.1 * log(f_difference[0]) + 0.05 * log(f_difference[1]) - 0.05 * log(f_difference[3]) - 0.1 * log(f_difference[4])) /
										(0.1429 * log(f_difference[0]) - 0.0714 * log(f_difference[1]) - 0.1429 * log(f_difference[2]) - 0.0714 * log(f_difference[3]) + 0.1429 * log(f_difference[4]));
									if (_isnan(delta))
									{
										delta = (0.5 * log(f_difference[1]) - 0.5 * log(f_difference[3])) / 2 /
											(0.5 * log(f_difference[1]) - log(f_difference[2]) + 0.5 * log(f_difference[3]));
									}
									if (_isnan(delta))
									{
										continue;
									}

									new_contour_points.push_back(cv::Point2f(x + delta * cos(sita), y + delta * sin(sita)));
								}

								//曲率滤波
								std::vector<float> curvature_vector;
								CalCurvatureFromEdgePoints(new_contour_points, curvature_vector);

								cv::RotatedRect new_rRect = fitEllipse(new_contour_points);    //fitEllipse只接受float和int类型

								//尝试通过拟合误差剔除坏点
						// 		ReduceBadEllipseFitPoints( new_contour_points,new_rRect.center.x,new_rRect.center.y,
						// 			new_rRect.size.width*0.5,new_rRect.size.height*0.5,new_rRect.angle*M_PI/180);
						// 		new_rRect =fitEllipse(new_contour_points);

								sub_pixel_center_x = new_rRect.center.x;
								sub_pixel_center_y = new_rRect.center.y;

								return true;
							}
							else
								if (subPixel_method == Gray_Moment)
								{
									///灰度矩方法，一维连续函数前3阶灰度矩
									std::vector<cv::Point2f> new_contour_points;

									for (int n = 0; n < contour_points.size(); n++)
									{
										int y = contour_points[n].y;
										int x = contour_points[n].x;

										if (x - 3 < 0 || x + 3 > image_mat.cols - 1 || y - 3 < 0 || y + 3 > image_mat.rows - 1)
										{
											continue;
										}

										if (x - center_x == 0)
										{
											continue;
										}

										float sita = atan2(y - center_y, x - center_x);
										float k = tan(sita);
										//float k = (y-center_y)/(x-center_x);   //梯度方向
										// 			double sita = gray_gradient_sita.at<double>(y,x);
										// 			double k = tan(sita);
										//double f_a,f_b,f_c,f_d;
										float abs_k = abs(k);

										float positive_gray[3];
										float negative_gray[3];

										int invt = 1;
										if (k > 0)
										{
											invt = -1;
										}
										if (abs_k < 1)
										{
											for (int i = 1; i < 4; i++)
											{
												int a = int(i * abs_k) / 1;
												float lamd2 = i * abs_k - a;
												float lamd1 = 1 - lamd2;

												positive_gray[i - 1] = lamd1 * image_mat.at<uchar>(y - a, x + invt * i) + lamd2 * image_mat.at<uchar>(y - a - 1, x + invt * i);
												negative_gray[i - 1] = lamd1 * image_mat.at<uchar>(y + a, x - invt * i) + lamd2 * image_mat.at<uchar>(y + a + 1, x - invt * i);
											}
										}
										else
										{
											for (int i = 1; i < 4; i++)
											{
												int a = int(i / abs_k) / 1;
												float lamd2 = i / abs_k - a;
												float lamd1 = 1 - lamd2;

												positive_gray[i - 1] = lamd1 * image_mat.at<uchar>(y - i, x + invt * a) + lamd2 * image_mat.at<uchar>(y - i, x + invt * (a + 1));
												negative_gray[i - 1] = lamd1 * image_mat.at<uchar>(y + i, x - invt * a) + lamd2 * image_mat.at<uchar>(y + i, x - invt * (a + 1));
											}
										}

										float gray_value[7];    //沿灰度梯度方向 ，向外
										if (y < center_y)
										{
											gray_value[0] = negative_gray[2];
											gray_value[1] = negative_gray[1];
											gray_value[2] = negative_gray[0];
											gray_value[3] = image_mat.at<uchar>(y, x);
											gray_value[4] = positive_gray[0];
											gray_value[5] = positive_gray[1];
											gray_value[6] = positive_gray[2];
										}
										else
										{
											gray_value[0] = positive_gray[2];
											gray_value[1] = positive_gray[1];
											gray_value[2] = positive_gray[0];
											gray_value[3] = image_mat.at<uchar>(y, x);
											gray_value[4] = negative_gray[0];
											gray_value[5] = negative_gray[1];
											gray_value[6] = negative_gray[2];
										}

										if (gray_value[0] > gray_value[6])
										{
											std::swap(gray_value[0], gray_value[6]);
											std::swap(gray_value[1], gray_value[5]);
											std::swap(gray_value[2], gray_value[4]);
											sita = sita + M_PI;
										}

										////计算灰度矩
										int N = 7;
										float m1 = 0, m2 = 0, m3 = 0;
										for (int i = 0; i < N; i++)
										{
											m1 += gray_value[i];
											m2 += gray_value[i] * gray_value[i];
											m3 += gray_value[i] * gray_value[i] * gray_value[i];
										}
										m1 /= N;
										m2 /= N;
										m3 /= N;

										float sigm = sqrt(m2 - m1 * m1);
										float s = (m3 + 2 * m1 * m1 * m1 - 3 * m1 * m2) / (sigm * sigm * sigm);

										float delta = N * 0.5 * s * sqrt(1.0 / (4 + s * s)) + (N + 1) * 0.5 - (N / 2 + 1);
										if (_isnan(delta))
										{
											continue;
										}

										new_contour_points.push_back(cv::Point2f(x + delta * cos(sita), y + delta * sin(sita)));
									}
									cv::RotatedRect new_rRect = fitEllipse(new_contour_points);    //fitEllipse只接受float和int类型

									//尝试通过拟合误差剔除坏点
									ReduceBadEllipseFitPoints(new_contour_points, new_rRect.center.x, new_rRect.center.y,
										new_rRect.size.width * 0.5, new_rRect.size.height * 0.5, new_rRect.angle * M_PI / 180);
									new_rRect = fitEllipse(new_contour_points);

									sub_pixel_center_x = new_rRect.center.x;
									sub_pixel_center_y = new_rRect.center.y;

									return true;
								}

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

float ImageDetectMethod::DICOTSU20140215(const cv::Mat& ori_image, int total_level /*=256*/)
{
	cv::Mat src;
	ori_image.convertTo(src, CV_16U);     //CV_16U 对应ushort

	int image_width = ori_image.cols;
	int image_height = ori_image.rows;

	int N = total_level;           //整个图分的等级
	//int pn[N] ={0};                //直方图的比例
	std::vector<int> pn;
	for (int i = 0; i < N; i++)
	{
		pn.push_back(0);
	}

	for (int i = 0; i < image_height; i++)
	{
		ushort* _src = src.ptr<ushort>(i);
		for (int j = 0; j < image_width; j++)
		{
			pn[int(_src[j])]++;
		}
	}

	float scale = 1. / (image_width * image_height);
	float mean = 0;
	for (int i = 0; i < N; i++)
	{
		mean += i * (float)pn[i];
	}
	mean *= scale;

	float q1 = 0, mean1 = 0;
	float max_sigma = 0, max_val = 0;

	for (int i = 0; i < N; i++)
	{
		float p_i, q2, mean2, sigma;

		p_i = pn[i] * scale;
		mean1 *= q1;
		q1 += p_i;
		q2 = 1. - q1;

		if (std::min(q1, q2) < FLT_EPSILON || std::max(q1, q2) > 1. - FLT_EPSILON)
			continue;

		mean1 = (mean1 + i * p_i) / q1;
		mean2 = (mean - q1 * mean1) / q2;

		sigma = q1 * q2 * (mean1 - mean2) * (mean1 - mean2);
		//sigma = q1*(mean1 - mean)*(mean1 - mean) + q2*(mean2 - mean)*(mean2 - mean);

		if (sigma > max_sigma)
		{
			max_sigma = sigma;
			max_val = i;
		}
	}

	return max_val;
}

bool ImageDetectMethod::CalCentriodBySubsetMat(const cv::Mat& subset_mat, float gray_threshold, float& sub_center_x, float& sub_center_y,
	MarkPointColorType color_type /*= BlackDownWhiteUp*/,
	SubPixelPosMethod subPixel_method /*= Gray_Centroid*/)
{
	sub_center_x = 0;
	sub_center_y = 0;
	float all_weight = 0;

	for (int i = 0; i < subset_mat.rows; i++)
	{
		const uchar* _subset_mat = subset_mat.ptr<uchar>(i);
		for (int j = 0; j < subset_mat.cols; j++)
		{
			float weight_value = CalWeightOfCentriod(float(_subset_mat[j]), gray_threshold, color_type, subPixel_method);
			sub_center_x += j * weight_value;
			sub_center_y += i * weight_value;
			all_weight += weight_value;
		}
	}

	if (all_weight == 0)
	{
		sub_center_x = 0;
		sub_center_y = 0;
	}
	else
	{
		sub_center_x /= all_weight;
		sub_center_y /= all_weight;
	}

	return true;
}

float ImageDetectMethod::CalWeightOfCentriod(float gray_value, float gray_threshold,
	MarkPointColorType color_type /*= BlackDownWhiteUp*/,
	SubPixelPosMethod subPixel_method /*= Gray_Centroid*/)
{
	//gray_threshold = 100;
	float weight_value;
	switch (color_type)
	{
	case BlackDownWhiteUp:
	{
		switch (subPixel_method)
		{
		case Binary_Centroid:
			if (gray_value - gray_threshold > 0)
			{
				weight_value = 1.;
			}
			else
				weight_value = 0;
			break;

		case Gray_Centroid:
			if (gray_value - gray_threshold > 0)
			{
				weight_value = gray_value /*- gray_threshold*/;
			}
			else
				weight_value = 0;

			//真实重心法
			//weight_value = gray_value;
			break;

		case Squared_Gray_Centroid:
			if (gray_value - gray_threshold > 0)
			{
				weight_value = (gray_value /*- gray_threshold*/) * (gray_value /*- gray_threshold*/);
			}
			else
				weight_value = 0;

			//真实灰度平方重心法5
			//weight_value = gray_value*gray_value;
			break;
		}
		break;
	}

	case WhiteDownBlackUp:
	{
		switch (subPixel_method)
		{
		case Binary_Centroid:
			if (gray_threshold - gray_value > 0)
			{
				weight_value = 1.;
			}
			else
				weight_value = 0;
			break;

		case Gray_Centroid:
			if (gray_threshold - gray_value > 0)
			{
				weight_value = gray_threshold - gray_value;
			}
			else
				weight_value = 0;
			break;

		case Squared_Gray_Centroid:
			if (gray_threshold - gray_value > 0)
			{
				weight_value = (gray_threshold - gray_value) * (gray_threshold - gray_value);
			}
			else
				weight_value = 0;
			break;
		}
		break;
	}
	}

	return weight_value;
}

MarkPointColorType ImageDetectMethod::JudgeTargetColorType(const cv::Mat& sub_mat, float center_x_insubmat, float center_y_insubmat,
	float ellipse_a, float ellipse_b, float angle_in_pi)
{
	float target_gray = 0;
	float back_ground_gray = 0;
	int target_pixel_num = 0;
	int back_ground_pixel_num = 0;
	for (int i = 0; i < sub_mat.rows; i++)
	{
		const uchar* _sub_mat_ptr = sub_mat.ptr<uchar>(i);
		for (int j = 0; j < sub_mat.cols; j++)
		{
			//坐标转换图像坐标系转椭圆坐标系
			float tr_x = (j - center_x_insubmat) * cos(angle_in_pi) + (i - center_y_insubmat) * sin(angle_in_pi);
			float tr_y = -(j - center_x_insubmat) * sin(angle_in_pi) + (i - center_y_insubmat) * cos(angle_in_pi);
			if (tr_x * tr_x / ellipse_a / ellipse_a + tr_y * tr_y / ellipse_b / ellipse_b < 1)
			{
				target_gray += _sub_mat_ptr[j];
				target_pixel_num++;
			}
			else
			{
				back_ground_gray += _sub_mat_ptr[j];
				back_ground_pixel_num++;
			}
		}
	}
	target_gray /= target_pixel_num;
	back_ground_gray /= back_ground_pixel_num;

	float gray_thresh = (target_gray + back_ground_gray) / 2;

	if (target_gray > back_ground_gray)
	{
		return BlackDownWhiteUp;
	}
	else
		return WhiteDownBlackUp;
}

float ImageDetectMethod::CalThresholdInSubsetMat(const cv::Mat& sub_mat, float center_x_insubmat, float center_y_insubmat,
	float ellipse_a, float ellipse_b, float angle_in_pi)
{
	//背景图和目标图灰度的均值
	float target_gray = 0;
	float back_ground_gray = 0;
	int target_pixel_num = 0;
	int back_ground_pixel_num = 0;
	for (int i = 0; i < sub_mat.rows; i++)
	{
		const uchar* _sub_mat_ptr = sub_mat.ptr<uchar>(i);
		for (int j = 0; j < sub_mat.cols; j++)
		{
			//坐标转换图像坐标系转椭圆坐标系
			float tr_x = (j - center_x_insubmat) * cos(angle_in_pi) + (i - center_y_insubmat) * sin(angle_in_pi);
			float tr_y = -(j - center_x_insubmat) * sin(angle_in_pi) + (i - center_y_insubmat) * cos(angle_in_pi);
			if (tr_x * tr_x / ellipse_a / ellipse_a + tr_y * tr_y / ellipse_b / ellipse_b < 1)
			{
				target_gray += _sub_mat_ptr[j];
				target_pixel_num++;
			}
			else
			{
				back_ground_gray += _sub_mat_ptr[j];
				back_ground_pixel_num++;
			}
		}
	}
	target_gray /= target_pixel_num;
	back_ground_gray /= back_ground_pixel_num;

	float gray_thresh = (target_gray + back_ground_gray) / 2;

	return gray_thresh;
}

float ImageDetectMethod::CalThresholdInSubsetMat2(const cv::Mat& image_mat, const std::vector<cv::Point>& contour_points)
{
	float threshold_value = 0;
	for (int i = 0; i < contour_points.size(); i++)
	{
		threshold_value += image_mat.at<uchar>(contour_points[i].y, contour_points[i].x);
	}
	return threshold_value / contour_points.size();
}
void ImageDetectMethod::CalGrayGradientBySobel(const cv::Mat& image_mat, int i, int j, float* gray_gradient, float* gradient_sita/*=NULL*/)
{
	//计算灰度梯度，这里通过sobel算子计算梯度幅值，和角度
	//sobel  x[-1, 0, 1          y[-1 -1 -1
	//         -1, 0, 1             0, 0, 0
	//         -1, 0, 1]            1, 1, 1]

	const uchar* _up_row_ptr = image_mat.ptr<uchar>(i - 1);
	const uchar* _mid_row_ptr = image_mat.ptr<uchar>(i);
	const uchar* _down_row_ptr = image_mat.ptr<uchar>(i + 1);
	float dx = (float(_up_row_ptr[j + 1] - _up_row_ptr[j - 1]) + 2 * float(_mid_row_ptr[j + 1] - _mid_row_ptr[j - 1]) + float(_down_row_ptr[j + 1] - _down_row_ptr[j - 1])) / 3;
	float dy = (float(_down_row_ptr[j - 1] - _up_row_ptr[j - 1]) + 2 * float(_down_row_ptr[j] - _up_row_ptr[j]) + float(_down_row_ptr[j + 1] - _up_row_ptr[j + 1])) / 3;

	*gray_gradient = sqrt(dx * dx + dy * dy);
	if (gradient_sita != NULL)
	{
		*gradient_sita = atan2(dy, dx);
	}
}

void ImageDetectMethod::CalCurvatureFromEdgePoints(const std::vector<cv::Point2f>& edge_points, std::vector<float>& curvature_vector)
{
	//<<光学测量中椭圆圆心定位算法研究>>  张虎， 2008
	//<<The Pre-Processing of Data Points for Curve Fitting in Reverse>>  2000
	//
	curvature_vector.clear();

	if (edge_points.size() < 3)
	{
		return;
	}

	float x1, y1, x2, y2, x3, y3, a, b, c, d, e, f, g, x0, y0, k;
	for (int i = 0; i < edge_points.size(); i++)
	{
		if (i == 0)
		{
			x1 = edge_points.at(edge_points.size() - 1).x;
			y1 = edge_points.at(edge_points.size() - 1).y;
			x2 = edge_points.at(i).x;
			y2 = edge_points.at(i).y;
			x3 = edge_points.at(i + 1).x;
			y3 = edge_points.at(i + 1).y;
		}
		else if (i == edge_points.size() - 1)
		{
			x1 = edge_points.at(i - 1).x;
			y1 = edge_points.at(i - 1).y;
			x2 = edge_points.at(i).x;
			y2 = edge_points.at(i).y;
			x3 = edge_points.at(0).x;
			y3 = edge_points.at(0).y;
		}
		else
		{
			x1 = edge_points.at(i - 1).x;
			y1 = edge_points.at(i - 1).y;
			x2 = edge_points.at(i).x;
			y2 = edge_points.at(i).y;
			x3 = edge_points.at(i + 1).x;
			y3 = edge_points.at(i + 1).y;
		}

		a = (x1 + x2) * (x2 - x1) * (y3 - y2);
		b = (x2 + x3) * (x3 - x2) * (y2 - y1);
		c = (y1 - y3) * (y2 - y1) * (y3 - y2);
		d = 2 * ((x2 - x1) * (y3 - y2) - (x3 - x2) * (y2 - y1));
		e = (y1 + y2) * (y2 - y1) * (x3 - x2);
		f = (y2 + y3) * (y3 - y2) * (x2 - x1);
		g = (x1 - x3) * (x2 - x1) * (x3 - x2);
		x0 = (a - b + c) / d;
		y0 = -(e - f + g) / d;

		k = 1. / sqrt((x0 - x2) * (x0 - x2) + (y0 - y2) * (y0 - y2));

		curvature_vector.push_back(k);
	}
}

void ImageDetectMethod::ReduceBadEllipseFitPoints(std::vector<cv::Point2f>& edge_points, float center_x, float center_y,
	float ellipse_a, float ellipse_b, float angle_in_pi)
{
	//尝试通过拟合误差剔除坏点
	std::vector<float> error_vector;
	for (int i = 0; i < edge_points.size(); i++)
	{
		float error = ErrorDROfEllipseFit(center_x, center_y, ellipse_a, ellipse_b, angle_in_pi,
			edge_points[i].x, edge_points[i].y);

		error_vector.push_back(error);

		if (error > 0.5)
		{
			std::vector<cv::Point2f>::iterator it = edge_points.begin() + i;
			edge_points.erase(it);

			i--;
		}
	}
}

void ImageDetectMethod::ShowContous(int image_width, int image_height,
	std::vector<std::vector<cv::Point>>& contours, QList<QList<float>>& ellipse_pars)
{
	cv::Mat show_contours_image = 255 * cv::Mat::ones(image_height, image_width, CV_8U);
	for (int i = 0; i < ellipse_pars.size(); i++)
	{
		drawContours(show_contours_image, contours, ellipse_pars[i][6], cv::Scalar(0));
	}

	//cv::namedWindow("canny", cv::WINDOW_FREERATIO);
	//imshow("canny", show_contours_image);
	////imwrite("E:\\canny.bmp",show_contours_image);
	//cv::waitKey(0);
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

double ImageDetectMethod::OTSUForCannyAdaptiveHighTh(cv::Mat ori_image, int total_level)
{
	cv::Mat src;
	ori_image.convertTo(src, CV_16U);     //CV_16U 对应ushort

	int image_width = ori_image.cols;
	int image_height = ori_image.rows;

	int N = total_level;           //整个图分的等级
	//int pn[N] ={0};                //直方图的比例
	std::vector<int> pn;
	for (int i = 0; i < N; i++)
	{
		pn.push_back(0);
	}

	for (int i = 0; i < image_height; i++)
	{
		ushort* _src = src.ptr<ushort>(i);
		for (int j = 0; j < image_width; j++)
		{
			pn[int(_src[j])]++;
		}
	}

	double scale = 1. / (image_width * image_height);
	double mean = 0;
	for (int i = 0; i < N; i++)
	{
		mean += i * (double)pn[i];
	}
	mean *= scale;

	double q1 = 0, mean1 = 0;
	double max_sigma = 0, max_val = 0;

	for (int i = 0; i < N; i++)
	{
		double p_i, q2, mean2, sigma;

		p_i = pn[i] * scale;
		mean1 *= q1;
		q1 += p_i;
		q2 = 1. - q1;

		if (std::min(q1, q2) < FLT_EPSILON || std::max(q1, q2) > 1. - FLT_EPSILON)
			continue;

		mean1 = (mean1 + i * p_i) / q1;
		mean2 = (mean - q1 * mean1) / q2;

		sigma = q1 * q2 * (mean1 - mean2) * (mean1 - mean2);
		//sigma = q1*(mean1 - mean)*(mean1 - mean) + q2*(mean2 - mean)*(mean2 - mean);

		if (i < N / 8)
		{
			continue;
		}
		if (sigma > max_sigma)
		{
			max_sigma = sigma;
			max_val = i;
		}
	}

	return max_val;
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
//*************************************
bool ImageDetectMethod::TestCircleCenterPosAccuracy(QString image_file_name, cv::Mat& code_point_mat, cv::Mat& uncode_point_mat,
	float ratio_k, float ratio_k1, float ratio_k2, float min_radius, float max_radius
	, float ellipse_error_pixel /*=0.5*/,
	MarkPointColorType color_type /*= BlackDownWhiteUp*/, CodePointBitesType code_bites_type /*=CodeBites15*/,
	DetectContoursMethod image_process_method /*= ADAPTIVE_THRESH_Method*/, SubPixelPosMethod subpixel_pos_method /*=Gray_Centroid*/)
{
	image_process_method = CANNY_Method;
	subpixel_pos_method = Gray_Moment;
	ImageDetectMethod::CodeAndUncodePointDetect(image_file_name, code_point_mat, uncode_point_mat,
		ratio_k, ratio_k1, ratio_k2,
		min_radius, max_radius, ellipse_error_pixel,
		color_type, code_bites_type,
		image_process_method, subpixel_pos_method);

	QFileInfo file_info(image_file_name);
	//WirteMatToFile(uncode_point_mat,file_info.completeBaseName());

	return true;
}

void ImageDetectMethod::WirteMatToFile(cv::Mat mat, QString file_name)
{
	QFile save_file("E:\\" + file_name + ".csv");
	if (save_file.open(QIODevice::WriteOnly | QIODevice::Text))
	{
		QTextStream stream(&save_file);

		if (mat.depth() == CV_64F)
		{
			for (int i = 0; i < mat.rows; i++)
			{
				for (int j = 0; j < mat.cols; j++)
				{
					stream << mat.at<double>(i, j) << ",";
				}
				stream << "\n";
			}
		}
		else if (mat.depth() == CV_8U)
		{
			for (int i = 0; i < mat.rows; i++)
			{
				for (int j = 0; j < mat.cols; j++)
				{
					stream << mat.at<uchar>(i, j) << ",";
				}
				stream << "\n";
			}
		}
		else if (mat.depth() == CV_32F)
		{
			for (int i = 0; i < mat.rows; i++)
			{
				for (int j = 0; j < mat.cols; j++)
				{
					stream << mat.at<float>(i, j) << ",";
				}
				stream << "\n";
			}
		}
		else if (mat.depth() == CV_16SC1)
		{
			for (int i = 0; i < mat.rows; i++)
			{
				for (int j = 0; j < mat.cols; j++)
				{
					stream << mat.at<short int>(i, j) << ",";
				}
				stream << "\n";
			}
		}
		else if (mat.depth() == CV_32S)
		{
			for (int i = 0; i < mat.rows; i++)
			{
				for (int j = 0; j < mat.cols; j++)
				{
					stream << mat.at<int>(i, j) << ",";
				}
				stream << "\n";
			}
		}
	}

	save_file.close();
}
//***********************************************************
bool ImageDetectMethod::FindChessGridOpencv(QString image_file_name, int h_num, int v_num, std::vector<cv::Point2f>& corners, std::vector<uchar>& sign_list, int& useful_corner_num)
{
	cv::Mat ori_image = cv::imread(image_file_name.toStdString(), 0);
	if (!ori_image.data)        // 判断图片调入是否成功
		return false;        // 调入图片失败则退出

	sign_list.clear();
	sign_list.resize(h_num * v_num);

	bool is_find = findChessboardCorners(ori_image, cv::Size(h_num, v_num), corners, cv::CALIB_CB_ADAPTIVE_THRESH + cv::CALIB_CB_NORMALIZE_IMAGE + cv::CALIB_CB_FAST_CHECK);

	if (is_find)
	{
		cornerSubPix(ori_image, corners, cv::Size(11, 11), cv::Size(-1, -1), cv::TermCriteria(cv::TermCriteria::COUNT + cv::TermCriteria::EPS, 30, 0.1));

		for (int i = 0; i < h_num * v_num; i++)
		{
			sign_list.at(i) = 1;
		}
		useful_corner_num = h_num * v_num;
		return true;
	}
	else
		return false;
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
							//FindSubPixelPosOfCircleCenter20140210(processed_image_mat, ellipse_pars[n][0], ellipse_pars[n][1], ellipse_pars[n][2],
							//	ellipse_pars[n][3], ellipse_pars[n][4], contours_copy[ellipse_pars[n][6]], sub_pixel_x, sub_pixel_y,
							//	&edge_contour, subpixel_pos_method, Uncertainty);
							contours_for_key[(i - 1) * h_num + j - 1] = contours_copy[ellipse_pars[n][6]];
							//corners[(i - 1) * h_num + j - 1] = cv::Point2f(ellipse_pars[n][0], ellipse_pars[n][1]);
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
							//FindSubPixelPosOfCircleCenter20140210(processed_image_mat, ellipse_pars[n][0], ellipse_pars[n][1], ellipse_pars[n][2],
							//	ellipse_pars[n][3], ellipse_pars[n][4], contours_copy[ellipse_pars[n][6]], sub_pixel_x, sub_pixel_y,
							//	&edge_contour, subpixel_pos_method, Uncertainty);
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

bool ImageDetectMethod::FindCodePointsForSelfCalibrationByHartly(QStringList image_file_name_list, std::vector<std::vector<cv::Point2f>>& code_points_list,
	float ratio_k /*= 2*/, float ratio_k1 /*= 2.4*/, float ratio_k2 /*= 4 */,
	float min_radius, float max_radius, float ellipse_error_pixel /*=0.5*/,
	MarkPointColorType color_type /*= BlackDownWhiteUp*/,
	CodePointBitesType code_bites_type /*=CodeBites15*/,
	DetectContoursMethod image_process_method /*= OTSU_Method*/,
	SubPixelPosMethod subpixel_pos_method /*=Gray_Centroid*/)
{
	if (image_file_name_list.size() < 4)
	{
		return false;
	}

	QList<cv::Mat> code_point_mat_list;
	for (int i = 0; i < image_file_name_list.size(); i++)
	{
		cv::Mat code_point_mat, uncode_point_mat;
		if (CodeAndUncodePointDetect(image_file_name_list[i], code_point_mat, uncode_point_mat,
			ratio_k, ratio_k1, ratio_k2,
			min_radius, max_radius, ellipse_error_pixel,
			color_type, code_bites_type,
			image_process_method, subpixel_pos_method))
		{
			code_point_mat_list.append(code_point_mat);
		}
	}

	code_points_list.clear();
	code_points_list.resize(4);

	/*vector<vector<cv::Point2f>> point_list(4);*/
	for (int i = 0; i < code_point_mat_list[0].rows; i++)
	{
		int match_index[4] = { -1,-1,-1,-1 };
		match_index[0] = i;
		for (int j = 1; j < 4; j++)
		{
			for (int k = 0; k < code_point_mat_list[j].rows; k++)
			{
				if (code_point_mat_list[0].at<float>(i, 0) == code_point_mat_list[j].at<float>(k, 0))
				{
					match_index[j] = k;
					break;
				}
			}
			if (match_index[j] == -1)
			{
				break;
			}
			else if (j == 3)
			{
				code_points_list[0].push_back(cv::Point2f(code_point_mat_list[0].at<float>(match_index[0], 1), code_point_mat_list[0].at<float>(match_index[0], 2)));
				code_points_list[1].push_back(cv::Point2f(code_point_mat_list[1].at<float>(match_index[1], 1), code_point_mat_list[1].at<float>(match_index[1], 2)));
				code_points_list[2].push_back(cv::Point2f(code_point_mat_list[2].at<float>(match_index[2], 1), code_point_mat_list[2].at<float>(match_index[2], 2)));
				code_points_list[3].push_back(cv::Point2f(code_point_mat_list[3].at<float>(match_index[3], 1), code_point_mat_list[3].at<float>(match_index[3], 2)));
			}
		}
	}

	if (code_points_list[0].size() < 9)
	{
		return false;
	}
	else
		return true;
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
				FindSubPixelPosOfCircleCenter20140210(processed_image_mat, ellipse_pars[i][0], ellipse_pars[i][1], ellipse_pars[i][2],
					ellipse_pars[i][3], ellipse_pars[i][4], contours[ellipse_pars[i][6]], sub_pixel_x, sub_pixel_y,
					&edge_contour, subpixel_pos_method, color_type);
				ellipse_pars[i][0] = sub_pixel_x;
				ellipse_pars[i][1] = sub_pixel_y;
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
				FindSubPixelPosOfCircleCenter20140210(processed_image_mat, ellipse_pars[i][0], ellipse_pars[i][1], ellipse_pars[i][2],
					ellipse_pars[i][3], ellipse_pars[i][4], contours[ellipse_pars[i][6]], sub_pixel_x, sub_pixel_y,
					&edge_contour, subpixel_pos_method, color_type);
				ellipse_pars[i][0] = sub_pixel_x;
				ellipse_pars[i][1] = sub_pixel_y;
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
			FindSubPixelPosOfCircleCenter20140210(processed_image_mat, ellipse_pars[i][0], ellipse_pars[i][1], ellipse_pars[i][2],
				ellipse_pars[i][3], ellipse_pars[i][4], contours[ellipse_pars[i][6]], sub_pixel_x, sub_pixel_y,
				&edge_contour, subpixel_pos_method, color_type);
			ellipse_pars[i][0] = sub_pixel_x;
			ellipse_pars[i][1] = sub_pixel_y;
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
			FindSubPixelPosOfCircleCenter20140210(processed_image_mat, ellipse_pars[i][0], ellipse_pars[i][1], ellipse_pars[i][2],
				ellipse_pars[i][3], ellipse_pars[i][4], contours[ellipse_pars[i][6]], sub_pixel_x, sub_pixel_y,
				&edge_contour, subpixel_pos_method, color_type);
			ellipse_pars[i][0] = sub_pixel_x;
			ellipse_pars[i][1] = sub_pixel_y;
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