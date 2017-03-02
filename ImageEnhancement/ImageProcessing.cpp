#include "stdafx.h"
#include "ImageProcessing.h"


CImageProcessing::CImageProcessing()
{
}


CImageProcessing::~CImageProcessing()
{
}


void CImageProcessing::retinex(cv::Mat src, cv::Mat im_dst, RetinexParams retinex_param)
{
	
	int width = src.cols;
	int height = src.rows;
	int size;
	int img_channel = src.channels();

	int          scale, row, col;
	int          i, j;
	int          pos;
	int          channel;
	uchar       *psrc = NULL;           
	float       *dst = nullptr;            
	float       *pdst = NULL;           
	float       *in, *out;
	int          channelsize;           
	float        weight;
	gauss_param  param;
	float        mean, var;
	float        mini, range, maxi;
	float        alpha;
	float        gain;
	float        offset;
	double       max_preview = 0.0;
	float        scales[6];

	size = width * height * src.channels();
	dst = (float*)GlobalAlloc(GMEM_FIXED | GMEM_ZEROINIT, size * sizeof(float));
	if (dst == nullptr)
		return;


	channelsize = (width * height);
	in = (float *)GlobalAlloc(GMEM_FIXED | GMEM_ZEROINIT, channelsize * sizeof(float));
	if (in == nullptr)
	{
		GlobalFree(dst);
		return; 
	}

	out = (float *)GlobalAlloc(GMEM_FIXED | GMEM_ZEROINIT, channelsize * sizeof(float));
	if (out == NULL)
	{
		GlobalFree(in);
		GlobalFree(dst);
		return;
	}

	scales_alloc(scales, retinex_param.scale, retinex_param.nscales, retinex_param.scales_mode);

	weight = 1. / retinex_param.nscales;

	pos = 0;
	for (channel = 0; channel < img_channel; channel++)
	{
		for (int i = 0; i < height; i++)
		{
			uchar* pS = src.ptr<uchar>(i);
			for (int j = 0; j < width; j++)
			{
				in[i * width + j] = pS[j * img_channel + channel] + 1.0;
			}
		}
		for (scale = 0; scale < retinex_param.nscales; scale++)
		{
			param_compute(&param, scales[scale]);
			
			for (row = 0; row < height; row++)
			{
				pos = row * width;
				gaussSmooth(in + pos, out + pos, width, 1, &param);
			}

			for (col = 0; col < width; col++)
			{
				pos = col;
				gaussSmooth(out + pos, in + pos, height, width, &param);
			}

			for (i = 0, pos = channel; i < channelsize; i++, pos += img_channel)
			{
				dst[pos] += weight * (log(src.data[pos] + 1.) - log(in[i]));
			}

		}
	}

	GlobalFree(in);
	GlobalFree(out);

	alpha = 128.;
	gain = 1.;
	offset = 0.;

	for (i = 0; i < size; i += img_channel)
	{
		float logl;

		psrc = src.data + i;
		pdst = dst + i;

		logl = log((float)psrc[0] + (float)psrc[1] + (float)psrc[2] + 3.);

		pdst[0] = gain * ((log(alpha * (psrc[0] + 1.)) - logl) * pdst[0]) + offset;
		pdst[1] = gain * ((log(alpha * (psrc[1] + 1.)) - logl) * pdst[1]) + offset;
		pdst[2] = gain * ((log(alpha * (psrc[2] + 1.)) - logl) * pdst[2]) + offset;
	}

	pdst = dst;

	compute_mean_var(pdst, &mean, &var, size, img_channel);
	mini = mean - 2*var;
	maxi = mean + 2*var;
	range = maxi - mini;

	if (!range)
		range = 1.0;

	for (i = 0; i < size; i += img_channel)
	{
		psrc = im_dst.data + i;
		pdst = dst + i;

		for (j = 0; j < 3; j++)
		{
			float c = 255 * (pdst[j] - mini) / range;

			if (c < 0)
			{
				c = 0;
			}
			else if (c > 255)
			{
				c = 255;
			}

			psrc[j] = (uchar)c;
		}
	}

	GlobalFree(dst);

}


void CImageProcessing::scales_alloc(float* scales, int scale, int nscales, int mode)
{
	if (nscales == 1)
	{
		scales[0] = (int)scale / 2;
	}
	else if (nscales == 2)
	{ 
		scales[0] = (int)scale / 2;
		scales[1] = (int)scale;
	}
	else
	{
		float size_step = (float)scale / (float)nscales;
		int   i;

		switch (mode)
		{
			case RETINEX_UNIFORM:
				for (i = 0; i < nscales; ++i)
					scales[i] = scale - (float)i * size_step;
				break;

			case RETINEX_LOW:
				size_step = (float)log(scale - 2.0) / (float)nscales;
				for (i = 0; i < nscales; ++i)
					scales[i] = pow(10, (i * size_step) / log(10));
				break;

			case RETINEX_HIGH:
				size_step = (float)log(scale - 2.0) / (float)nscales;
				for (i = 0; i < nscales; ++i)
					scales[i] = scale - pow(10, (i * size_step) / log(10));
				break;

			default:
				break;
		}
	}
}


void CImageProcessing::param_compute(gauss_param* params, float sigma)
{
	float q, q2, q3;
	q = 0;
	if (sigma >= 2.5)
	{
		q = 0.98711 * sigma - 0.96330;
	}
	else if ((sigma >= 0.5) && (sigma < 2.5))
	{
		q = 3.97156 - 4.14554 * (float)sqrt((double)1 - 0.26891 * sigma);
	}
	else
	{
		q = 0.1147705018520355224609375;
	}

	q2 = q * q;
	q3 = q * q2;
	params->b[0] = (1.57825 + (2.44413*q) + (1.4281 *q2) + (0.422205*q3));
	params->b[1] = ((2.44413*q) + (2.85619*q2) + (1.26661 *q3));
	params->b[2] = (-((1.4281*q2) + (1.26661 *q3)));
	params->b[3] = ((0.422205*q3));
	params->B = 1.0 - ((params->b[1] + params->b[2] + params->b[3]) / params->b[0]);
	params->sigma = sigma;
}

void CImageProcessing::gaussSmooth(float* in, float* out, int size, int stride, gauss_param *param)
{
	int i, n, bufsize;
	float *w1, *w2;

	/* 正向迭代 */
	bufsize = size + 3;
	size -= 1;
	w1 = (float *)GlobalAlloc(GMEM_ZEROINIT | GMEM_FIXED, bufsize * sizeof(float));
	w2 = (float *)GlobalAlloc(GMEM_FIXED | GMEM_ZEROINIT, bufsize * sizeof(float));
	w1[0] = in[0];
	w1[1] = in[0];
	w1[2] = in[0];
	for (i = 0, n = 3; i <= size; i++, n++)
	{
		w1[n] = (float)(param->B * in[i * stride] +
				((param->b[1] * w1[n - 1] +
				  param->b[2] * w1[n - 2] +
				  param->b[3] * w1[n - 3] ) / param->b[0]));
	}

	/* 逆向迭代 */
	w2[size + 1] = w1[size + 3];
	w2[size + 2] = w1[size + 3];
	w2[size + 3] = w1[size + 3];
	for (i = size, n = i; i >= 0; i--, n--)
	{
		w2[n] = out[i * stride] = (float)(param->B * w1[n] +
				((param->b[1] * w2[n + 1] +
				  param->b[2] * w2[n + 2] +
				  param->b[3] * w2[n + 3] ) / param->b[0]));
	}

	GlobalFree(w1);
	GlobalFree(w2);
}


void CImageProcessing::compute_mean_var(float *src, float *mean, float *var, int size, int nchannel)
{
	float vsquared;
	int i, j;
	float *psrc;

	vsquared = 0;
	*mean = 0;
	for (i = 0; i < size; i += nchannel)
	{
		psrc = src + i;
		for (j = 0; j < nchannel; j++)
		{
			*mean += psrc[j];
			vsquared += psrc[j] * psrc[j];
		}
	}

	*mean /= (float)size; /* mean */
	vsquared /= (float)size; /* mean (x^2) */
	*var = (vsquared - (*mean * *mean));
	*var = sqrt(*var); /* var */
}