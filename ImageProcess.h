#pragma once

#include <atlimage.h>
#include <map>
#include <cmath>
#include <algorithm>

#include "OpenCLUtils.h"

#define MAX_THREAD 8
#define MAX_SPAN 15

#define calcX(x) x % distW
#define calcY(y) y / distW
#define index(x, y, w) y * w + x
#define srcPoint(x, y) srcData + pit*(y)+(x)*cnt
#define distPoint(x, y) distData + dPit*(y)+(x)*dCnt
#define getColor(r, g, b, p) r = *(p + 2), g = *(p + 1), b = *p

using namespace std;

typedef unsigned long long ull;

// 图元
struct ImageFragment {
	CImage *src;
	CImage *dist;
	double start;
	double end;
	int maxSpan; // 为模板中心到边缘的距离
	int index; // 区域编号
};

struct NoiseParam {
	int type = 0; // 0: 无噪声, 1: 高斯, 2: 椒盐
	double avg = 0, std = 0, factor = 0;
};

struct FilterParam {
	int type = 0; // 0: 无滤波, 1: 中值, 2: 均值, 3: 高斯, 4: 维纳, 5: 双边
	double std = 0;
};

// 处理参数
struct ProcessParam {
	double angle, xScale, yScale;
	int algo, threadCount;
	double startTime;
	bool dft = false;

	NoiseParam noise;
	FilterParam filter;
};

// 点
struct Point {
	double x, y;
	Point() : x(0), y(0) {}
	Point(double x, double y) : x(x), y(y) {}
};

static class ImageProcess {
public:
	static const double PI;
	static const double PI2;
	static const double E;

	static const double E2PI;

	static const double FOURIER_FACTOR;

	static const double SALT_NOISE_FACTOR;

	static const double GAUS_NOISE_AVG;
	static const double GAUS_NOISE_STD;

	// 中值滤波核大小（需为奇数）
	static const int MEDIAN_KERNEL_SIZE;
	// 均值滤波核大小
	static const int MEAN_KERNEL_SIZE;
	// 高斯滤波核大小
	static const int GAUS_KERNEL_SIZE;
	// 维纳滤波核大小
	static const int WIENER_KERNEL_SIZE;
	// 双边滤波核大小
	static const int BILATERAL_KERNEL_SIZE;

	static void initialize();

	static CImage* createImage(CImage* src);

	static ImageFragment* createImageFragment(CImage *src, CImage *dist);
	static ImageFragment* createImageFragment(CImage *src, CImage *dist, int count, int index = 0);
	static ImageFragment* createImageFragment(CImage *src, CImage *dist, 
		double start, double end, int index, int maxSpan);

	static void process(ImageFragment *frag, bool msg = true);
	static void processNoise(ImageFragment *frag);
	static void processFilter(ImageFragment *frag);

	static void setParam(ProcessParam param);

	static void copyImage(CImage* src, CImage* dist);
	static void copyImage(ImageFragment *frag);

	static void translateImage(ImageFragment * frag, 
		double xScale, double yScale, double angle, int algo = 0);

	static void translateImage(ImageFragment * frag);

	static void dftTranslate(ImageFragment * frag);

	//static void clDftTranslate(ImageFragment * frag);

	static void saltNoise(ImageFragment * frag, double factor);
	static void saltNoise(ImageFragment * frag);

	static void gaussianNoise(ImageFragment * frag, double avg, double std);
	static void gaussianNoise(ImageFragment * frag);

	static void meanFilter(ImageFragment * frag);

	static void medianFilter(ImageFragment * frag);

	static void medianFilterCL(ImageFragment * frag);

	static void gaussianFilter(ImageFragment * frag, double std);
	static void gaussianFilter(ImageFragment * frag);

	static void wienerFilter(ImageFragment * frag);

	static void bilateralFiter(ImageFragment * frag, double std);
	static void bilateralFiter(ImageFragment * frag);

	/*
	static void rotateImage(ImageFragment *frag, double angle);
	static void rotateImage(ImageFragment *frag);

	static void scaleImage(ImageFragment *frag,
		double xScale, double yScale, int mode = 0);
	static void scaleImage(ImageFragment *frag);
	*/

private:
	// static map<int, double> bcCache;
	static ProcessParam param;

	static CSize calcTranslatedSize(CImage *src);
	/*
	static double _biCubicPoly(double x);
	static void _biCubicScale(ImageFragment *frag, double xScale, double yScale);
	static void _normalScale(ImageFragment *frag, double xScale, double yScale);
	*/
	
	static double _biCubicPoly(double x);

	static void _biCubicTranslate(ull distS, ull distE, int distW, int distH,
		int distW2, int distH2, float fCosa, float fSina, double mx, double my,
		byte *srcData, int pit, int cnt, int srcW, int srcH, byte * distData, int dPit, int dCnt);

	static void _biCubicTranslateCL(byte *srcData, byte *distData, int distW, int distH,
		int srcW, int srcH, int pit, int dPit, double angle, double xScale, double yScale);
	
	static void _normalTranslate(ull distS, ull distE, int distW, int distH,
		int distW2, int distH2, float fCosa, float fSina, double mx, double my, 
		byte *srcData, int pit, int cnt, int srcW, int srcH, byte * distData, int dPit, int dCnt);

	//static void _2DArrayScale(ImageFragment *frag, double xScale, double yScale);

	static double _generateBoxMuller(double mean, double stddev);

	static void _generateGaussianKernel(double t[3][3], double std);
	static void _generateGaussianKernelForBin(double t[3][3], double std);
	static void _generateColorKernel(double c[256], double std);

	static size_t roundUp(int groupSize, int globalSize);
};
