//#include "stdafx.h"
//#include "ImageProcessing.h"
//
//
//
//float *dImage = NULL;   //original image
//float *dResult = NULL;   //temp array for iterations
//size_t pitch;
//int radius;
//
//
//
//void compute_mean_var(float *src, float *mean, float *var, int size, int nchannel)
//{
//	float vsquared;
//	int i, j;
//	float *psrc;
//
//	vsquared = 0;
//	*mean = 0;
//	for (i = 0; i < size; i += nchannel)
//	{
//		psrc = src + i;
//		for (j = 0; j < nchannel; j++)
//		{
//			*mean += psrc[j];
//			vsquared += psrc[j] * psrc[j];
//		}
//	}
//
//	*mean /= (float)size; /* mean */
//	vsquared /= (float)size; /* mean (x^2) */
//	*var = (vsquared - (*mean * *mean));
//	*var = sqrt(*var); /* var */
//}
//
//void gaussSmooth(float* in, float* out, int size, int stride, gauss_param *param)
//{
//	int i, n, bufsize;
//	float *w1, *w2;
//
//	/* 正向迭代 */
//	bufsize = size + 3;
//	size -= 1;
//	w1 = (float *)GlobalAlloc(GMEM_ZEROINIT | GMEM_FIXED, bufsize * sizeof(float));
//	w2 = (float *)GlobalAlloc(GMEM_FIXED | GMEM_ZEROINIT, bufsize * sizeof(float));
//	w1[0] = in[0];
//	w1[1] = in[0];
//	w1[2] = in[0];
//	for (i = 0, n = 3; i <= size; i++, n++)
//	{
//		w1[n] = (float)(param->B * in[i * stride] +
//				((param->b[1] * w1[n - 1] +
//				  param->b[2] * w1[n - 2] +
//				  param->b[3] * w1[n - 3] ) / param->b[0]));
//	}
//
//	/* 逆向迭代 */
//	w2[size + 1] = w1[size + 3];
//	w2[size + 2] = w1[size + 3];
//	w2[size + 3] = w1[size + 3];
//	for (i = size, n = i; i >= 0; i--, n--)
//	{
//		w2[n] = out[i * stride] = (float)(param->B * w1[n] +
//				((param->b[1] * w2[n + 1] +
//				  param->b[2] * w2[n + 2] +
//				  param->b[3] * w2[n + 3] ) / param->b[0]));
//	}
//
//	GlobalFree(w1);
//	GlobalFree(w2);
//}
//
//void param_compute(gauss_param* params, float sigma)
//{
//	float q, q2, q3;
//	q = 0;
//	if (sigma >= 2.5)
//	{
//		q = 0.98711 * sigma - 0.96330;
//	}
//	else if ((sigma >= 0.5) && (sigma < 2.5))
//	{
//		q = 3.97156 - 4.14554 * (float)sqrt((double)1 - 0.26891 * sigma);
//	}
//	else
//	{
//		q = 0.1147705018520355224609375;
//	}
//
//	q2 = q * q;
//	q3 = q * q2;
//	params->b[0] = (1.57825 + (2.44413*q) + (1.4281 *q2) + (0.422205*q3));
//	params->b[1] = ((2.44413*q) + (2.85619*q2) + (1.26661 *q3));
//	params->b[2] = (-((1.4281*q2) + (1.26661 *q3)));
//	params->b[3] = ((0.422205*q3));
//	params->B = 1.0 - ((params->b[1] + params->b[2] + params->b[3]) / params->b[0]);
//	params->sigma = sigma;
//}
//
//void scales_alloc(float* scales, int scale, int nscales, int mode)
//{
//	if (nscales == 1)
//	{
//		scales[0] = (int)scale / 2;
//	}
//	else if (nscales == 2)
//	{ 
//		scales[0] = (int)scale / 2;
//		scales[1] = (int)scale;
//	}
//	else
//	{
//		float size_step = (float)scale / (float)nscales;
//		int   i;
//
//		switch (mode)
//		{
//			case RETINEX_UNIFORM:
//				for (i = 0; i < nscales; ++i)
//					scales[i] = scale - (float)i * size_step;
//				break;
//
//			case RETINEX_LOW:
//				size_step = (float)log(scale - 2.0) / (float)nscales;
//				for (i = 0; i < nscales; ++i)
//					scales[i] = pow(10, (i * size_step) / log(10));
//				break;
//
//			case RETINEX_HIGH:
//				size_step = (float)log(scale - 2.0) / (float)nscales;
//				for (i = 0; i < nscales; ++i)
//					scales[i] = scale - pow(10, (i * size_step) / log(10));
//				break;
//
//			default:
//				break;
//		}
//	}
//}
//
//void retinex(cv::Mat src, cv::Mat im_dst, RetinexParams retinex_param)
//{
//	int width = src.cols;
//	int height = src.rows;
//	int size;
//	int img_channel = src.channels();
//
//	int          scale, row, col;
//	int          i, j;
//	int          pos;
//	int          channel;
//	uchar       *psrc = NULL;           
//	float       *dst = nullptr;            
//	float       *pdst = NULL;           
//	float       *in, *out;
//	int          channelsize;           
//	float        weight;
//	gauss_param  param;
//	float        mean, var;
//	float        mini, range, maxi;
//	float        alpha;
//	float        gain;
//	float        offset;
//	double       max_preview = 0.0;
//	float        scales[6];
//
//	size = width * height * src.channels();
//	dst = (float*)GlobalAlloc(GMEM_FIXED | GMEM_ZEROINIT, size * sizeof(float));
//	if (dst == nullptr)
//		return;
//
//
//	channelsize = (width * height);
//	in = (float *)GlobalAlloc(GMEM_FIXED | GMEM_ZEROINIT, channelsize * sizeof(float));
//	if (in == nullptr)
//	{
//		GlobalFree(dst);
//		return; 
//	}
//
//	out = (float *)GlobalAlloc(GMEM_FIXED | GMEM_ZEROINIT, channelsize * sizeof(float));
//	if (out == NULL)
//	{
//		GlobalFree(in);
//		GlobalFree(dst);
//		return;
//	}
//
//	scales_alloc(scales, retinex_param.scale, retinex_param.nscales, retinex_param.scales_mode);
//
//	weight = 1. / retinex_param.nscales;
//
//	pos = 0;
//	for (channel = 0; channel < img_channel; channel++)
//	{
//		for (int i = 0; i < height; i++)
//		{
//			uchar* pS = src.ptr<uchar>(i);
//			for (int j = 0; j < width; j++)
//			{
//				in[i * width + j] = pS[j * img_channel + channel] + 1.0;
//			}
//		}
//		for (scale = 0; scale < retinex_param.nscales; scale++)
//		{
//			param_compute(&param, scales[scale]);
//			// 待优化
//			for (row = 0; row < height; row++)
//			{
//				pos = row * width;
//				gaussSmooth(in + pos, out + pos, width, 1, &param);
//			}
//
//			for (col = 0; col < width; col++)
//			{
//				pos = col;
//				gaussSmooth(out + pos, in + pos, height, width, &param);
//			}
//
//			for (i = 0, pos = channel; i < channelsize; i++, pos += img_channel)
//			{
//				dst[pos] += weight * (log(src.data[pos] + 1.) - log(in[i]));
//			}
//
//		}
//	}
//
//	GlobalFree(in);
//	GlobalFree(out);
//
//	alpha = 128.;
//	gain = 1.;
//	offset = 0.;
//
//	for (i = 0; i < size; i += img_channel)
//	{
//		float logl;
//
//		psrc = src.data + i;
//		pdst = dst + i;
//
//		logl = log((float)psrc[0] + (float)psrc[1] + (float)psrc[2] + 3.);
//
//		pdst[0] = gain * ((log(alpha * (psrc[0] + 1.)) - logl) * pdst[0]) + offset;
//		pdst[1] = gain * ((log(alpha * (psrc[1] + 1.)) - logl) * pdst[1]) + offset;
//		pdst[2] = gain * ((log(alpha * (psrc[2] + 1.)) - logl) * pdst[2]) + offset;
//	}
//
//	pdst = dst;
//
//	compute_mean_var(pdst, &mean, &var, size, img_channel);
//	mini = mean - 2*var;
//	maxi = mean + 2*var;
//	range = maxi - mini;
//
//	if (!range)
//		range = 1.0;
//
//	for (i = 0; i < size; i += img_channel)
//	{
//		psrc = im_dst.data + i;
//		pdst = dst + i;
//
//		for (j = 0; j < 3; j++)
//		{
//			float c = 255 * (pdst[j] - mini) / range;
//
//			if (c < 0)
//			{
//				c = 0;
//			}
//			else if (c > 255)
//			{
//				c = 255;
//			}
//
//			psrc[j] = (uchar)c;
//		}
//	}
//
//	GlobalFree(dst);
//
//}
//
//
//
//void getT(float* dark, float* t, int width, int height, int r)
//{
//	int d = 2 * r + 1;
//	float *left2right = (float *)::GlobalAlloc(GMEM_ZEROINIT | GMEM_FIXED, sizeof(float)*width);
//	float *right2left = (float *)::GlobalAlloc(GMEM_ZEROINIT | GMEM_FIXED, sizeof(float)*width);
//	float *high2low = (float *)::GlobalAlloc(GMEM_ZEROINIT | GMEM_FIXED, sizeof(float)*height);
//	float *low2high = (float *)::GlobalAlloc(GMEM_ZEROINIT | GMEM_FIXED, sizeof(float)*height);
//
//	// 行方向
//	for (int i = 0; i < height; i++)
//	{
//		for (int j = 0; j < width; j += d)
//		{
//			for (int k = 0; k < d && j + k < width; k++)
//			{
//				int left = j + k;
//				int right = j + d - k - 1;
//				if (right >= width)
//				{
//					right = j + width%d - k - 1;
//				}
//				if (k != 0)
//				{
//					left2right[left] = dark[i*width + left] < left2right[left - 1] ? dark[i*width + left] : left2right[left - 1];
//					right2left[right] = dark[i*width + right] < right2left[right + 1] ? dark[i*width + right] : right2left[right + 1];
//				}
//				else
//				{
//					left2right[left] = dark[i*width + left];
//					right2left[right] = dark[i*width + right];
//				}
//			}
//		}
//		for (int j = 0; j < width; j++)
//		{
//			if (j < r)
//			{
//				t[i*width + j] = left2right[j + r];
//			}
//			else if (width - j - 1 < r)
//			{
//				t[i*width + j] = right2left[j - r];
//			}
//			else
//			{
//				t[i*width + j] = left2right[j + r] < right2left[j - r] ? left2right[j + r] : right2left[j - r];
//			}
//		}
//	}
//	// 列方向
//	for (int i = 0; i < width; i++)
//	{
//		for (int j = 0; j < height; j += d)
//		{
//			for (int k = 0; k < d && j + k < height; k++)
//			{
//				int low = j + k;
//				int high = j + d - k - 1;
//				if (high >= height)
//				{
//					high = j + height%d - k - 1;
//				}
//				if (k != 0)
//				{
//					low2high[low] = t[low*width + i] < low2high[low - 1] ? t[low*width + i] : low2high[low - 1];
//					high2low[high] = t[high*width + i] < high2low[high + 1] ? t[high*width + i] : high2low[high + 1];
//				}
//				else
//				{
//					low2high[low] = t[low*width + i];
//					high2low[high] = t[high*width + i];
//				}
//			}
//		}
//		for (int j = 0; j < height; j++)
//		{
//			if (j < r)
//			{
//				t[j*width + i] = low2high[j + r];
//			}
//			else if (height - j - 1 < r)
//			{
//				t[j*width + i] = high2low[j - r];
//			}
//			else
//			{
//				t[j*width + i] = low2high[j + r] < high2low[j - r] ? low2high[j + r] : high2low[j - r];
//			}
//		}
//	}
//	::GlobalFree(right2left);
//	::GlobalFree(left2right);
//	::GlobalFree(high2low);
//	::GlobalFree(low2high);
//}
//
//void boxFilter(float* src, int width, int height, int radius)
//{
//	float *buff = (float*)::GlobalAlloc(GMEM_FIXED | GMEM_ZEROINIT, sizeof(float)*width*height);
//	// 列方向
//	for (int i = 0; i < height; i++)
//	{
//		if (i != 0)
//			for (int j = 0; j < width; j++)
//				buff[i*width + j] = src[i*width + j] + buff[(i - 1)*width + j];
//		else
//			for (int j = 0; j < width; j++)
//				buff[i*width + j] = src[i*width + j];
//	}
//	for (int i = 0; i < radius + 1; i++)
//	{
//		for (int j = 0; j < width; j++)
//		{
//			src[i*width + j] = buff[(radius + i)*width + j];
//		}
//	}
//	for (int i = radius + 1; i < height - radius - 1; i++)
//	{
//		for (int j = 0; j < width; j++)
//		{
//			src[i*width + j] = buff[(radius + i)*width + j] - buff[(i - radius - 1)*width + j];
//		}
//	}
//	for (int i = height - radius - 1; i < height; i++)
//	{
//		for (int j = 0; j < width; j++)
//		{
//			src[i*width + j] = buff[(height - 1)*width + j] - buff[(i - radius - 1)*width + j];
//		}
//	}
//	// 行方向
//	for (int i = 0; i < width; i++)
//	{
//		if (i != 0)
//			for (int j = 0; j < height; j++)
//				buff[j*width + i] = src[j*width + i] + buff[j*width + i - 1];
//		else
//			for (int j = 0; j < height; j++)
//				buff[j*width + i] = src[j*width + i];
//	}
//	for (int i = 0; i < radius + 1; i++)
//	{
//		for (int j = 0; j < height; j++)
//		{
//			src[j*width + i] = buff[j*width + (i + radius)];
//		}
//	}
//	for (int i = radius + 1; i < width - radius - 1; i++)
//	{
//		for (int j = 0; j < height; j++)
//		{
//			src[j*width + i] = buff[j*width + (i + radius)] - buff[j*width + (i - radius - 1)];
//		}
//	}
//	for (int i = width - radius - 1; i < width; i++)
//	{
//		for (int j = 0; j < height; j++)
//		{
//			src[j*width + i] = buff[j*width + width - 1] - buff[j*width + (i - radius - 1)];
//		}
//	}
//
//	int w = radius * 2 + 1;
//	for (int i = 0; i < radius; i++)
//	{
//		for (int j = i; j < width - i; j++)
//		{
//			if (j < radius)
//			{
//				src[i*width + j] /= ((radius + 1 + i)*(radius + 1 + j));
//				src[(height - i - 1)*width + j] /= ((radius + 1 + i)*(radius + 1 + j));
//			}
//			else if (width - j - 1 < radius)
//			{
//				src[i*width + j] /= ((radius + 1 + i)*(radius + 1 + width - j - 1));
//				src[(height - i - 1)*width + j] /= ((radius + 1 + i)*(radius + 1 + width - j - 1));
//			}
//			else
//			{
//				src[i*width + j] /= (radius + 1 + i)*w;
//				src[(height - i - 1)*width + j] /= (radius + 1 + i)*w;
//			}
//		}
//	}
//	for (int i = 0; i < radius; i++)
//	{
//		for (int j = i + 1; j < height - i - 1; j++)
//		{
//			if (j < radius)
//			{
//				src[j*width + i] /= ((radius + 1 + i)*(radius + 1 + j));
//				src[j*width + (width - i - 1)] /= ((radius + 1 + i)*(radius + 1 + j));
//			}
//			else if (height - j - 1 < radius)
//			{
//				src[j*width + i] /= ((radius + 1 + i)*(radius + 1 + height - j - 1));
//				src[j*width + (width - i - 1)] /= ((radius + 1 + i)*(radius + 1 + height - j - 1));
//			}
//			else
//			{
//				src[j*width + i] /= (radius + 1 + i)*w;
//				src[j*width + (width - 1 - i)] /= (radius + i + 1)*w;
//			}
//		}
//	}
//	for (int i = radius; i < height - radius; i++)
//	{
//		for (int j = radius; j < width - radius; j++)
//		{
//			src[i*width + j] /= w*w;
//		}
//	}
//	::GlobalFree(buff);
//}
//
//void guidedFilter(cv::Mat guide_image, float* src, int radius, double eps)
//{
//	int size = guide_image.cols * guide_image.rows;
//	int width = guide_image.cols, height = guide_image.rows;
//	float *temp_IP = (float*)::GlobalAlloc(GMEM_FIXED | GMEM_ZEROINIT, sizeof(float)*size);
//	float *temp_II = (float*)::GlobalAlloc(GMEM_FIXED | GMEM_ZEROINIT, sizeof(float)*size);
//	float *temp_P = (float*)::GlobalAlloc(GMEM_FIXED | GMEM_ZEROINIT, sizeof(float)*size);
//	float *temp_I = (float*)::GlobalAlloc(GMEM_FIXED | GMEM_ZEROINIT, sizeof(float)*size);
//	float *cov_IP = (float*)::GlobalAlloc(GMEM_FIXED | GMEM_ZEROINIT, sizeof(float)*size);
//	float *var_II = (float*)::GlobalAlloc(GMEM_FIXED | GMEM_ZEROINIT, sizeof(float)*size);
//	float *a = (float*)::GlobalAlloc(GMEM_FIXED | GMEM_ZEROINIT, sizeof(float)*size);
//	float *b = (float*)::GlobalAlloc(GMEM_FIXED | GMEM_ZEROINIT, sizeof(float)*size);
//
//	for (int i = 0; i < size; i++)
//	{
//		temp_P[i] = src[i];
//	}
//	for (int i = 0; i < size; i++)
//	{
//		temp_I[i] = guide_image.data[i];
//	}
//	for (int i = 0; i < size; ++i)
//	{
//		temp_IP[i] = temp_I[i] * temp_P[i];
//	}
//	for (int i = 0; i < size; ++i)
//	{
//		temp_II[i] = temp_I[i] * temp_I[i];
//	}
//	boxFilter(temp_I, width, height, radius);
//	boxFilter(temp_P, width, height, radius);
//	boxFilter(temp_IP, width, height, radius);
//	boxFilter(temp_II, width, height, radius);
//	for (int i = 0; i < size; ++i)
//	{
//		var_II[i] = temp_II[i] - temp_I[i] * temp_I[i];
//	}
//	for (int i = 0; i < size; ++i)
//	{
//		cov_IP[i] = temp_IP[i] - temp_I[i] * temp_P[i];
//	}
//	for (int i = 0; i < size; ++i)
//	{
//		a[i] = cov_IP[i] / (var_II[i] + eps);
//	}
//	for (int i = 0; i < size; ++i)
//	{
//		b[i] = temp_P[i] - a[i] * temp_I[i];
//	}
//	boxFilter(a, width, height, radius);
//	boxFilter(b, width, height, radius);
//	for (int i = 0; i < size; ++i)
//	{
//		guide_image.data[i] = a[i] * guide_image.data[i] + b[i];
//	}
//	::GlobalFree(temp_I);
//	::GlobalFree(temp_P);
//	::GlobalFree(temp_II);
//	::GlobalFree(temp_IP);
//	::GlobalFree(cov_IP);
//	::GlobalFree(var_II);
//	::GlobalFree(a);
//	::GlobalFree(b);
//	
//}
//
//void guideFilter_color(float* t, cv::Mat src, int radius, float eps)
//{
//	int  size = src.cols * src.rows;
//	cv::Mat  bgr[3];
//
//	cv::split(src, bgr);
//	for (int i = 0; i < 3; ++i)
//	{
//		guidedFilter(bgr[i], t, radius, eps);
//	}
//	for (int i = 0; i < size; i++)
//	{
//		t[i] = 0.299 * bgr[2].data[i] + 0.587 * bgr[1].data[i] + 0.114 * bgr[0].data[i];
//	}
//
//}
//
//void ImageProcessByDark(cv::Mat src, cv::Mat dst, GuideParames para)
//{
//	int width, height, size;
//	uchar       *psrc = nullptr;
//	float       *pdst = nullptr;
//	float       *dark = nullptr;//暗通道图
//	float       *t = nullptr;   //透射率图
//	int         i, j;
//	int         img_channel;
//	double sum = 0.;
//	int pointNum = 0;
//	float A;
//	double pix;
//	float stretch_p[256], stretch_p1[256], stretch_num[256];
//
//
//
//	width = src.cols;
//	height = src.rows;
//	size = width * height;
//	img_channel = src.channels();
//	assert(img_channel == 3);
//	dark = (float*)GlobalAlloc(GMEM_FIXED | GMEM_ZEROINIT, size * sizeof(float));
//	if (dark == nullptr)
//		return;
//
//	//获取暗原色
//	for (i = 0; i < size; ++i)
//	{
//		psrc = src.data + i*img_channel;
//		pdst = dark + i;
//
//		uchar b = psrc[0];
//		uchar g = psrc[1];
//		uchar r = psrc[2];
//		*pdst = ((b < g) ? b : g) < r ? ((b < g) ? b : g) : r;
//	}
//
//	
//
//	/*memset(stretch_p, 0, sizeof(stretch_p));
//	memset(stretch_p1, 0, sizeof(stretch_p1));
//	memset(stretch_num, 0, sizeof(stretch_num));
//	for (i = 0; i < size; ++i)
//	{
//		float pixel0 = dark[i];
//		int pixel = (int)pixel0;
//		stretch_num[pixel]++;
//	}
//	for (i = 0; i < 256; ++i)
//	{
//		stretch_p[i] = stretch_num[i] / size;
//	}
//	for (i = 0; i < 256; ++i)
//	{
//		for (j = 0; j <= i; ++j)
//		{
//			stretch_p1[i] += stretch_p[j];
//			if (stretch_p1[i] > 0.999)
//			{
//				pix = (double)i;
//				i = 256;
//				break;
//			}
//		}
//	}
//	for (i = 0; i < size; ++i)
//	{
//		double temp = dark[i];
//		if (temp > pix)
//		{
//			psrc = src.data + i*img_channel;
//			pointNum++;
//			sum += psrc[0];
//			sum += psrc[1];
//			sum += psrc[2];
//		}
//	}
//	A = sum / (3 * pointNum);
//	if (A > 220.0)
//	{*/
//		A = 250.0;
//	//}
//	
//	t = (float*)::GlobalAlloc(GMEM_FIXED | GMEM_ZEROINIT, size * sizeof(float));
//	if (t == nullptr)
//	{
//		::GlobalFree(dark);
//		return;
//	}
//	getT(dark, t, width, height, para.min_radius);
//	::GlobalFree(dark);
//
//	guideFilter_color(t, src, para.guide_radius, 0.0004);
//	for (int i = 0; i < size; ++i)
//	{
//		t[i] = 1 - 0.9 * (t[i] / A);
//		if (t[i] < 0)
//			t[i] = 0;
//	}
//	uchar* ppsrc;
//	for (i = 0; i < size; ++i)
//	{
//		ppsrc = src.data + i * img_channel;
//		psrc = dst.data + i*img_channel;
//		pdst = t + i;
//		if (*pdst < 0.1) *pdst = 0.1;
//		psrc[0] = (uchar)((ppsrc[0] - A*(1 - *pdst)) / *pdst);
//		psrc[1] = (uchar)((ppsrc[1] - A*(1 - *pdst)) / *pdst);
//		psrc[2] = (uchar)((ppsrc[2] - A*(1 - *pdst)) / *pdst);
//	}
//	::GlobalFree(t);
//}
//
//void getTmp(float* t, int width, int height, int r, float* tmp){
//	
//	int i = 0;
//	for (; i < r; ++i){
//		for (int j = 0; j < width; ++j){
//			tmp[i * (width + 2 * r) + (j + r)] = t[(r - i - 1) * width + j];
//		}
//	}
//	for (; i < height + r; ++i){
//		for (int j = 0; j < width; ++j){
//			tmp[i * (width + 2 * r) + (j + r)] = t[(i - r)*width + j];
//		}
//	}
//	for (; i < height + 2 * r; ++i){
//		for (int j = 0; j < width; ++j){
//			tmp[i * (width + 2 * r) + (j + r)] = t[(2*height - 1 - i + r)*width + j];
//		}
//	}
//	for (i = 0; i < height + 2 * r; ++i){
//		for (int j = 0; j < r; ++j){
//			tmp[i * (width + 2 * r) + (r - j - 1)] = tmp[i * (width + 2 * r) + (r + j)];
//		}
//	}
//	for (i = 0; i < height + 2 * r; ++i){
//		for (int j = 0; j < r; ++j){
//			tmp[i * (width + 2 * r) + (r + width + j)] = tmp[i * (width + 2 * r) + (r + width - j - 1)];
//		}
//	}
//}
//
//void BilateralFilter(float* t, int width, int height){
//
//	int d = 25;
//	double sigma_color = 40.0;
//	double sigma_space = 20;
//
//	int  maxk, radius;
//
//	if (sigma_color <= 0)
//		sigma_color = 1;
//	if (sigma_space <= 0)
//		sigma_space = 1;
//
//	double gauss_color_coeff = -0.5 / (sigma_color*sigma_color);
//	double gauss_space_coeff = -0.5 / (sigma_space*sigma_space);
//
//	if (d <= 0)
//		radius = cvRound(sigma_space*1.5);
//	else
//		radius = d / 2;
//	radius = MAX(radius, 1);
//	d = radius * 2 + 1;
//	float* tmp = (float*)GlobalAlloc(GMEM_FIXED | GMEM_ZEROINIT, (width + 2 * radius)*(height + 2 * radius)* sizeof(float));
//	getTmp(t, width, height, radius, tmp);
//	float* row_tmp = (float*)GlobalAlloc(GMEM_ZEROINIT | GMEM_FIXED, width*sizeof(float));
//	float* col_tmp = (float*)GlobalAlloc(GMEM_FIXED | GMEM_ZEROINIT, height*sizeof(float));
//
//
//	//颜色、空间权重
//	vector<float> _color_weight(256);
//	vector<float> _space_weight(d);
//	float* color_weight = &_color_weight[0];
//	float* space_weight = &_space_weight[0];
//	
//
//	int i;
//	for (i = 0; i < 256; i++)
//		color_weight[i] = (float)exp(i*i*gauss_color_coeff);
//
//	for (i = -radius; i <= radius; i++)
//	{
//		space_weight[i + radius] = (float)exp(i*i*gauss_space_coeff);
//	}
//
//	//水平方向双边滤波
//	for (int h = 0; h < height + 2*radius; h++)
//	{
//		float* psrc = tmp + h*(width + 2 * radius);
//		float* pdst = row_tmp;
//
//		for (int w = 0; w < width; w++)
//		{
//			float sum = 0, wsum = 0;
//			int val0 = psrc[w + radius];
//			int k = 0;
//			for (; k < d; k++)
//			{
//				int val = psrc[w + k];
//				float w = space_weight[k] * color_weight[abs(val - val0)];
//				sum += val*w;
//				wsum += w;
//			}
//			pdst[w] = (float)(sum / wsum);
//		}
//		for (int w = 0; w < width; ++w){
//			psrc[w + radius] = pdst[w];
//		}
//	}
//	//转置
//	float* tmp_transpose = (float*)GlobalAlloc(GMEM_ZEROINIT | GMEM_FIXED, (width + 2 * radius)*(height + 2 * radius)*sizeof(float));
//	for (int h = 0; h < width + 2 * radius; ++h){
//		for (int w = 0; w < height + 2 * radius; ++w){
//			tmp_transpose[h*(height + 2 * radius) + w] = tmp[w*(width + 2 * radius) + h];
//		}
//	}
//	GlobalFree(tmp);
//	//垂直方向双边滤波
//	for (int w = 0; w < width; w++)
//	{
//		float* psrc = tmp_transpose + (w + radius)*(height + 2 * radius);
//		float* pdst = col_tmp;
//		
//		for (int h = 0; h < height; h++)
//		{
//			float sum = 0, wsum = 0;
//			int val0 = psrc[h + radius];
//			int k = 0;
//			for (; k < d; k++)
//			{
//				int val = psrc[h + k];
//				float w = space_weight[k] * color_weight[abs(val - val0)];
//				sum += val*w;
//				wsum += w;
//			}
//			pdst[h] = (float)(sum / wsum);
//		}
//		for (int h = 0; h < height; ++h){
//			psrc[h + radius] = pdst[h];
//		}
//	}
//
//	for (int row = 0; row < height; ++row){
//		for (int col = 0; col < width; col++){
//			t[row*width + col] = tmp_transpose[(col + radius)*(height + 2 * radius) + (row + radius)];
//		}
//	}
//	GlobalFree(tmp_transpose);
//}
//
//void ImageProcessByDarkWithBilateralFilter(cv::Mat src, cv::Mat dst, GuideParames para){
//
//	int width, height, size;
//	uchar       *psrc = nullptr;
//	float       *pdst = nullptr;
//	float       *dark = nullptr;//暗通道图
//	float       *t = nullptr;   //透射率图
//	int         i, j;
//	int         img_channel;
//	double sum = 0.;
//	int pointNum = 0;
//	float A;
//	double pix;
//	float stretch_p[256], stretch_p1[256], stretch_num[256];
//
//
//
//	width = src.cols;
//	height = src.rows;
//	size = width * height;
//	img_channel = src.channels();
//	assert(img_channel == 3);
//	dark = (float*)GlobalAlloc(GMEM_FIXED | GMEM_ZEROINIT, size * sizeof(float));
//	if (dark == nullptr)
//		return;
//
//	//获取暗原色
//	for (i = 0; i < size; ++i)
//	{
//		psrc = src.data + i*img_channel;
//		pdst = dark + i;
//
//		uchar b = psrc[0];
//		uchar g = psrc[1];
//		uchar r = psrc[2];
//		*pdst = ((b < g) ? b : g) < r ? ((b < g) ? b : g) : r;
//	}
//	/*memset(stretch_p, 0, sizeof(stretch_p));
//	memset(stretch_p1, 0, sizeof(stretch_p1));
//	memset(stretch_num, 0, sizeof(stretch_num));
//	for (i = 0; i < size; ++i)
//	{
//		float pixel0 = dark[i];
//		int pixel = (int)pixel0;
//		stretch_num[pixel]++;
//	}
//	for (i = 0; i < 256; ++i)
//	{
//		stretch_p[i] = stretch_num[i] / size;
//	}
//	for (i = 0; i < 256; ++i)
//	{
//		for (j = 0; j <= i; ++j)
//		{
//			stretch_p1[i] += stretch_p[j];
//			if (stretch_p1[i] > 0.999)
//			{
//				pix = (double)i;
//				i = 256;
//				break;
//			}
//		}
//	}
//	for (i = 0; i < size; ++i)
//	{
//		double temp = dark[i];
//		if (temp > pix)
//		{
//			psrc = src.data + i*img_channel;
//			pointNum++;
//			sum += psrc[0];
//			sum += psrc[1];
//			sum += psrc[2];
//		}
//	}
//	A = sum / (3 * pointNum);
//	if (A > 220.0)
//	{*/
//		A = 250.0;
//	//}
//	t = (float*)::GlobalAlloc(GMEM_FIXED | GMEM_ZEROINIT, size * sizeof(float));
//	if (t == nullptr)
//	{
//		::GlobalFree(dark);
//		return;
//	}
//	getT(dark, t, width, height, para.min_radius);
//	::GlobalFree(dark);
//
//
//	BilateralFilter(t, width, height);
//	for (int i = 0; i < size; ++i)
//	{
//		t[i] = 1 - 0.9 * (t[i] / A);
//		if (t[i] < 0)
//			t[i] = 0;
//	}
//	uchar* ppsrc;
//	for (i = 0; i < size; ++i)
//	{
//		ppsrc = src.data + i * img_channel;
//		psrc = dst.data + i*img_channel;
//		pdst = t + i;
//		if (*pdst < 0.1) *pdst = 0.1;
//		psrc[0] = (uchar)((ppsrc[0] - A*(1 - *pdst)) / *pdst);
//		psrc[1] = (uchar)((ppsrc[1] - A*(1 - *pdst)) / *pdst);
//		psrc[2] = (uchar)((ppsrc[2] - A*(1 - *pdst)) / *pdst);
//	}
//	::GlobalFree(t);
//}
#include "stdafx.h"