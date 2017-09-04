#ifndef __Import_h__
#define __Import_h__
#include <cuda_runtime.h> 
#include "device_launch_parameters.h"

#include <opencv2/core/core.hpp>  
#include <opencv2/highgui/highgui.hpp>
#include <opencv/cv.h>
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <Windows.h>
#include <iostream>
using namespace std;

#define RETINEX_UNIFORM 0
#define RETINEX_LOW     1
#define RETINEX_HIGH    2
typedef struct
{
	int     scale;
	int     nscales;
	int     scales_mode;
} RetinexParams;

typedef struct
{
	int min_radius;
	int guide_radius;
}GuideParames;
#endif // __Import_h__