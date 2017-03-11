#pragma once

class CImageProcessing
{
public:
	CImageProcessing();
	~CImageProcessing();

private:
	typedef struct
	{
		double  sigma;
		double B;
		double b[4];
	} gauss_param;


public:
	/*
	 * 
	 */
	static void ImageProcessByDark(cv::Mat src, cv::Mat dst, GuideParames para);
	static void getT(float* dark, float* t, int width, int height, int r);
	static void guideFilter_color(float* T, cv::Mat src, int radius, float eps);
	static void guidedFilter(cv::Mat guide_image, float* src, int radius, double eps);
	static void boxFilter(float* src, int width, int height, int radius);
	/*
	 * 以下为retinex处理相关函数
	 */
	static void retinex(cv::Mat src, cv::Mat im_dst, RetinexParams retinex_param);
	static void scales_alloc(float* scales, int scale, int nscales, int mode);
	static void param_compute(gauss_param* params, float sigma);
	static void gaussSmooth(float* in, float* out, int size, int stride, gauss_param *param);
	static void compute_mean_var(float *src, float *mean, float *var, int size, int nchannel);
};

