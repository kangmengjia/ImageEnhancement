//#ifndef __ImageProcessing_h__
//#define __ImageProcessing_h__
//
//#include "stdafx.h"
//using namespace cv;
//
//
//
//typedef struct
//{
//	double  sigma;
//	double B;
//	double b[4];
//} gauss_param;
//
//void ImageProcessByDark(cv::Mat src, cv::Mat dst, GuideParames para);
//void ImageProcessByDarkWithBilateralFilter(cv::Mat src, cv::Mat dst, GuideParames para);
//void getTmp(float* t, int width, int height, int r, float* tmp);
//void BilateralFilter(float* t, int width, int height);
//void getT(float* dark, float* t, int width, int height, int r);
//void guideFilter_color(float* T, cv::Mat src, int radius, float eps);
//void guidedFilter(cv::Mat guide_image, float* src, int radius, double eps);
//void boxFilter(float* src, int width, int height, int radius);
///*
// * 以下为retinex处理相关函数
// */
//void retinex(cv::Mat src, cv::Mat im_dst, RetinexParams retinex_param);
//void scales_alloc(float* scales, int scale, int nscales, int mode);
//void param_compute(gauss_param* params, float sigma);
//void gaussSmooth(float* in, float* out, int size, int stride, gauss_param *param);
//void compute_mean_var(float *src, float *mean, float *var, int size, int nchannel);
//
//#endif // __ImageProcessing_h__
