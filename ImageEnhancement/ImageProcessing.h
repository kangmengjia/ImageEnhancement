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

	static void retinex(cv::Mat src, cv::Mat im_dst, RetinexParams retinex_param);
	static void scales_alloc(float* scales, int scale, int nscales, int mode);
	static void param_compute(gauss_param* params, float sigma);
	static void gaussSmooth(float* in, float* out, int size, int stride, gauss_param *param);
	static void compute_mean_var(float *src, float *mean, float *var, int size, int nchannel);
};

