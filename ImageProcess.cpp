#include "stdafx.h"
#include "Debug.h"
#include "ImageProcess.h"
#include "OpenCVMFCDlg.h"

#include <math.h>

const double ImageProcess::PI = 3.14159265358979;
const double ImageProcess::PI2 = 2*3.14159265358979;
const double ImageProcess::E = 2.718281828459;
const double ImageProcess::E2PI = pow(ImageProcess::E, 2 * ImageProcess::PI);

const double ImageProcess::FOURIER_FACTOR = 14;

const double ImageProcess::SALT_NOISE_FACTOR = 0.05;

const double ImageProcess::GAUS_NOISE_AVG = 50;
const double ImageProcess::GAUS_NOISE_STD = 100;

const int ImageProcess::MEDIAN_KERNEL_SIZE = 3;
const int ImageProcess::MEAN_KERNEL_SIZE = 3;
const int ImageProcess::GAUS_KERNEL_SIZE = 3;
const int ImageProcess::WIENER_KERNEL_SIZE = 3;
const int ImageProcess::BILATERAL_KERNEL_SIZE = 3;

ProcessParam ImageProcess::param;

void ImageProcess::process(ImageFragment *frag) {
	//LOG("Processing: " << frag->index); 
	translateImage(frag); 
	processNoise(frag); processFilter(frag);
	if (param.dft) dftTranslate(frag);
	::PostMessage(AfxGetMainWnd()->GetSafeHwnd(), WM_TRANSLATE, 1, NULL);
	delete frag;
}

void ImageProcess::processNoise(ImageFragment *frag) {
	switch (param.noise.type) {
	case 1: gaussianNoise(frag); break;
	case 2: saltNoise(frag); break;
	}
}

void ImageProcess::processFilter(ImageFragment *frag) {
	switch (param.filter.type) {
	case 1: medianFilter(frag); break;
	case 2: meanFilter(frag); break;
	case 3: gaussianFilter(frag); break;
	case 4: wienerFilter(frag); break;
	case 5: bilateralFiter(frag); break;
	}
}

void ImageProcess::setParam(ProcessParam param) {
	ImageProcess::param = param;
}

CImage* ImageProcess::createImage(CImage* src) {
	if (src != NULL) {
		src->Destroy(); delete src;
	}
	return new CImage();
}

ImageFragment* ImageProcess::createImageFragment(CImage *src, CImage *dist, 
	double start, double end, int index, int maxSpan) {
	ImageFragment *frag = new ImageFragment();

	CSize size = calcTranslatedSize(src);
	//LOG("Size: " << size.cx << ", " << size.cy);

	if (!dist->IsNull()) dist->Destroy();
	dist->Create(size.cx, size.cy, src->GetBPP(), 0);

	frag->src = src; frag->dist = dist;
	frag->start = start; frag->end = end;
	frag->index = index; frag->maxSpan = maxSpan;
	return frag;
}

ImageFragment* ImageProcess::createImageFragment(CImage *src, CImage *dist) {
	return createImageFragment(src, dist, 1);
}

ImageFragment* ImageProcess::createImageFragment(CImage *src, CImage *dist, int count, int index) {
	double rate = 1.0 / count, start = index * rate, end = start + rate;
	//LOG("rate: " << rate << ", start: " << start << ", end: " << end);
	return createImageFragment(src, dist, start, end, index, MAX_SPAN);
}

CSize ImageProcess::calcTranslatedSize(CImage *src) {
	int srcW = src->GetWidth(), srcH = src->GetHeight();
	double h = srcH * param.yScale, w = srcW * param.xScale;
	Point minP, maxP, tmpP;
	Point p1(-w/2, h/2), p2(w/2, h/2), p3(w/2, -h/2), p4(-w/2, -h/2);
	double d = sqrt(p1.x*p1.x + p1.y*p1.y);
	double a, r = param.angle * PI / 180;
	
	a = atan(p1.y / p1.x);
	minP = maxP = tmpP = Point(-d * cos(a + r), d * sin(a + r));

	a = atan(p2.y / p2.x);
	tmpP = Point(d * cos(a + r), d * sin(a + r));
	minP.x = min(minP.x, tmpP.x), minP.y = min(minP.y, tmpP.y);
	maxP.x = max(maxP.x, tmpP.x), maxP.y = max(maxP.y, tmpP.y);

	a = atan(p3.y / p3.x);
	tmpP = Point(d * cos(a + r), -d * sin(a + r));
	minP.x = min(minP.x, tmpP.x), minP.y = min(minP.y, tmpP.y);
	maxP.x = max(maxP.x, tmpP.x), maxP.y = max(maxP.y, tmpP.y);

	a = atan(p4.y / p4.x);
	tmpP = Point(-d * cos(a + r), -d * sin(a + r));
	minP.x = min(minP.x, tmpP.x), minP.y = min(minP.y, tmpP.y);
	maxP.x = max(maxP.x, tmpP.x), maxP.y = max(maxP.y, tmpP.y);

	return CSize(maxP.x - minP.x, maxP.y - minP.y);
}

void ImageProcess::copyImage(CImage* src, CImage* dist) {
	int srcH = src->GetHeight(), srcW = src->GetWidth();

	if (!dist->IsNull()) dist->Destroy();
	dist->Create(srcH, srcH, src->GetBPP(), 0);

	int MaxColors = src->GetMaxColorTableEntries();
	RGBQUAD* ColorTab;
	ColorTab = new RGBQUAD[MaxColors];

	CDC *pDCsrc, *pDCdrc;
	
	if (src->IsIndexed()) {
		src->GetColorTable(0, MaxColors, ColorTab);
		dist->SetColorTable(0, MaxColors, ColorTab);
	}

	pDCsrc = CDC::FromHandle(src->GetDC());
	pDCdrc = CDC::FromHandle(dist->GetDC());
	pDCdrc->BitBlt(0, 0, src->GetWidth(), src->GetHeight(), pDCsrc, 0, 0, SRCCOPY);
	src->ReleaseDC();
	dist->ReleaseDC();
	delete ColorTab;

}

void ImageProcess::copyImage(ImageFragment *frag) {
	CImage *src = frag->src, *dist = frag->dist;
	if (src == NULL || dist == NULL) return;
	if (src->IsNull()) return;
	copyImage(src, dist);
	/*
	int srcW = src->GetWidth(), srcH = src->GetHeight();
	int distW = srcW, distH = srcH;

	double start = frag->start, end = frag->end;

	int pit = src->GetPitch(), cnt = src->GetBPP() / 8;
	ull srcCnt = srcW * srcH, distCnt = distW * distH;
	ull srcS = srcCnt * start, srcE = srcCnt * end;

	//if (!dist->IsNull()) dist->Destroy();
	if (!dist->IsNull() && (dist->GetWidth() != distW || dist->GetHeight() != distH))
		dist->Destroy();
	if (dist->IsNull()) dist->Create(distW, distH, src->GetBPP(), 0);

	int dPit = dist->GetPitch(), dCnt = dist->GetBPP() / 8;

	//目标图像参数
	byte* srcData = (byte*)src->GetBits();
	byte* distData = (byte*)dist->GetBits();
	*/
	/*
	for (int i = 0; i < srcH; i++)
		memcpy(distData + i * dPit, srcData + i * pit, abs(pit));
	*/
	/*

	byte* p;

	int fr, fg, fb;

	//复制图像数据
	for (int i = 0; i < srcCnt; ++i) {
		int x = calcX(i), y = calcY(i);
		p = srcPoint(x, y);
		getColor(fr, fg, fb, p);

		*(distPoint(x, y) + 2) = fr;
		*(distPoint(x, y) + 1) = fg;
		*(distPoint(x, y)) = fb;
	}
	*/
}

void ImageProcess::translateImage(ImageFragment *frag, 
	double xScale, double yScale, double angle, int algo) {

	CImage *src = frag->src, *dist = frag->dist;
	if (src == NULL || dist == NULL) return;
	if (src->IsNull() || dist->IsNull()) return;

	int srcH = src->GetHeight(), srcW = src->GetWidth();
	int distH = dist->GetHeight(), distW = dist->GetWidth();
	
	int distH2 = srcH * yScale, distW2 = srcW * xScale; // 实际图片所占大小
	
	int pit = src->GetPitch(), cnt = src->GetBPP() / 8;
	int dPit = dist->GetPitch(), dCnt = dist->GetBPP() / 8;

	double start = frag->start, end = frag->end;

	ull srcCnt = srcW * srcH, distCnt = distW * distH;
	ull srcS = srcCnt * start, srcE = srcCnt * end;
	ull distS = distCnt * start, distE = distCnt * end;

	/*
	LOG("from " << start << " to " << end);
	LOG("srcS = " << srcS << ", srcE = " << srcE);
	LOG("distS = " << distS << ", distE = " << distE);
	LOG("srcCnt = " << srcCnt << ", distCnt = " << distCnt);
	LOG("distW = " << distW << ", distH = " << distH);
	LOG("distW2 = " << distW2 << ", distH2 = " << distH2);
	*/

	byte* srcData = (byte*)src->GetBits();
	byte* distData = (byte*)dist->GetBits();

	float fCosa = cos(angle * PI / 180), fSina = sin(angle * PI / 180);

	//LOG("fCosa = " << fCosa << ", fSina = " << fSina);

	double mx = ((double)distW2 - 1) / ((double)srcW - 1);
	double my = ((double)distH2 - 1) / ((double)srcH - 1);

	switch (algo) {
	case 0:
		_biCubicTranslate(distS, distE, distW, distH, distW2, distH2, fCosa, fSina,
			mx, my, srcData, pit, cnt, srcW, srcH, distData, dPit, dCnt);
		break;
	case 1:
		_normalTranslate(distS, distE, distW, distH, distW2, distH2, fCosa, fSina, 
			mx, my, srcData, pit, cnt, srcW, srcH, distData, dPit, dCnt);
		break;
	}
}

void ImageProcess::translateImage(ImageFragment * frag) {
	translateImage(frag, param.xScale, param.yScale, param.angle, param.algo);
}

void ImageProcess::dftTranslate(ImageFragment * frag) {
	CImage* dist = frag->dist;
	if (dist == NULL || dist->IsNull()) return;

	int distH = dist->GetHeight(), distW = dist->GetWidth();
	int dPit = dist->GetPitch(), dCnt = dist->GetBPP() / 8;

	double start = frag->start, end = frag->end;

	ull distCnt = distW * distH;
	ull distS = distCnt * start, distE = distCnt * end;

	byte* distData = (byte*)dist->GetBits();

	byte *p; // 临时点

	int u, v, x, y;
	double t;

	//int fr, fb, fg; // 最终颜色
	int tr, tb, tg; // 中间运算
	int grey;

	double real, imag, mag;

	for (ull i = distS; i < distE; ++i) {
		u = calcX(i), v = calcY(i);

		real = imag = 0;

		for (ull j = 0; j < distCnt; ++j) {
			x = calcX(j), y = calcY(j);
			p = distPoint(x, y);

			getColor(tr, tg, tb, p);

			t = cos(2 * PI * ((u*x*1.0) / distW + (v*y*1.0) / distH));
			grey = (tr * 30 + tg * 59 + tb * 11 + 50) / 100;

			if ((x + y) & 1) grey = -grey;

			real += grey * cos(t); imag -= grey * sin(t);
		}
		mag = sqrt(real * real + imag * imag);
		mag = FOURIER_FACTOR * log(mag + 1);
		mag = max(0, min(mag, 255));

		*(distPoint(u, v) + 2) = mag;
		*(distPoint(u, v) + 1) = mag;
		*(distPoint(u, v)) = mag;
	}
}

void ImageProcess::gaussianNoise(ImageFragment * frag, double avg, double std) {
	CImage* dist = frag->dist;
	if (dist == NULL || dist->IsNull()) return;

	int distH = dist->GetHeight(), distW = dist->GetWidth();
	int dPit = dist->GetPitch(), dCnt = dist->GetBPP() / 8;

	double start = frag->start, end = frag->end;

	ull distCnt = distW * distH;
	ull distS = distCnt * start, distE = distCnt * end;

	byte* distData = (byte*)dist->GetBits();

	byte *p; // 临时点

	int x, y;
	int fr, fb, fg; // 最终颜色

	double tmp;

	for (ull i = distS; i < distE; ++i) {
		x = calcX(i), y = calcY(i);
		p = distPoint(x, y);

		getColor(fr, fg, fb, p);

		tmp = fr + _generateBoxMuller(avg, std);
		fr = max(min(tmp, 255), 0);
		tmp = fg + _generateBoxMuller(avg, std);
		fg = max(min(tmp, 255), 0);
		tmp = fb + _generateBoxMuller(avg, std);
		fb = max(min(tmp, 255), 0);

		*(distPoint(x, y) + 2) = fr;
		*(distPoint(x, y) + 1) = fg;
		*(distPoint(x, y)) = fb;
	}
}

void ImageProcess::gaussianNoise(ImageFragment * frag) {
	gaussianNoise(frag, param.noise.avg, param.noise.std);
}

void ImageProcess::meanFilter(ImageFragment * frag) {

	CImage* dist = frag->dist;
	if (dist == NULL || dist->IsNull()) return;

	const int SIZE = MEAN_KERNEL_SIZE, N = SIZE * SIZE;
	const int D = SIZE >> 1, M = N >> 1;

	int distH = dist->GetHeight(), distW = dist->GetWidth();
	int dPit = dist->GetPitch(), dCnt = dist->GetBPP() / 8;

	double start = frag->start, end = frag->end;

	ull distCnt = distW * distH;
	ull distS = distCnt * start, distE = distCnt * end;

	byte* distData = (byte*)dist->GetBits();

	byte *p; // 临时点

	int x, y;
	int fr, fg, fb;
	int tr, tg, tb;

	// 要比较的像素值
	byte r[N], g[N], b[N];
	for (ull i = distS; i < distE; ++i) {
		x = calcX(i), y = calcY(i);
		p = distPoint(x, y);

		int index = 0, di = i - distS;

		// 处理边界
		if (x < D || x >= distW - D || y < D || y >= distH - D) 
			continue;

		fr = fg = fb = 0;
		// 获取周围 SIZE*SIZE 的颜色
		for (int dx = -D; dx <= D; ++dx)
			for (int dy = -D; dy <= D; ++dy) {
				p = distPoint(x + dx, y + dy);
				getColor(tr, tg, tb, p);
				fr += tr; fg += tg; fb += tb;
				if (index++ >= N) break;
			}

		// 求平均并赋值
		*(distPoint(x, y) + 2) = fr / N;
		*(distPoint(x, y) + 1) = fg / N;
		*(distPoint(x, y)) = fb / N;
	}
}

void ImageProcess::medianFilter(ImageFragment * frag) {
	CImage* dist = frag->dist;
	if (dist == NULL || dist->IsNull()) return;

	const int SIZE = MEDIAN_KERNEL_SIZE, N = SIZE * SIZE;
	const int D = SIZE >> 1, M = N >> 1;

	int distH = dist->GetHeight(), distW = dist->GetWidth();
	int dPit = dist->GetPitch(), dCnt = dist->GetBPP() / 8;

	double start = frag->start, end = frag->end;

	ull distCnt = distW * distH;
	ull distS = distCnt * start, distE = distCnt * end;
	ull length = distE - distS;

	byte* distData = (byte*)dist->GetBits();

	byte* fr = new byte[length];
	byte* fg = new byte[length];
	byte* fb = new byte[length];

	byte *p; // 临时点

	int x, y;

	// 要比较的像素值
	byte r[N], g[N], b[N];
	for (ull i = distS; i < distE; ++i) {
		x = calcX(i), y = calcY(i);
		p = distPoint(x, y);

		int index = 0, di = i - distS;

		// 处理边界
		if (x < D || x >= distW - D || y < D || y >= distH - D) {
			getColor(fr[di], fg[di], fb[di], p);
			continue;
		}

		// 获取周围 SIZE*SIZE 的颜色
		for (int dx = -D; dx <= D; ++dx)
			for (int dy = -D; dy <= D; ++dy) {
				p = distPoint(x + dx, y + dy);
				getColor(r[index], g[index], b[index], p);
				if (index++ >= N) break;
			}

		// 排序，求出中值
		sort(r, r + N); sort(g, g + N); sort(b, b + N);

		fr[di] = r[M], fg[di] = g[M], fb[di] = b[M];
	}

	for (ull i = 0; i < length; ++i) {
		x = calcX(i), y = calcY(i);
		*(distPoint(x, y) + 2) = fr[i];
		*(distPoint(x, y) + 1) = fg[i];
		*(distPoint(x, y)) = fb[i];
	}

	delete[] fr;
	delete[] fg;
	delete[] fb;
}

void ImageProcess::gaussianFilter(ImageFragment * frag, double std) {

	CImage* dist = frag->dist;
	if (dist == NULL || dist->IsNull()) return;

	const int SIZE = GAUS_KERNEL_SIZE, N = SIZE * SIZE;
	const int D = SIZE >> 1, M = N >> 1;

	int distH = dist->GetHeight(), distW = dist->GetWidth();
	int dPit = dist->GetPitch(), dCnt = dist->GetBPP() / 8;

	double start = frag->start, end = frag->end;

	ull distCnt = distW * distH;
	ull distS = distCnt * start, distE = distCnt * end;

	byte* distData = (byte*)dist->GetBits();

	byte *p; // 临时点

	int x, y;
	double fr, fg, fb;
	int tr, tg, tb;

	double mult;
	double m[SIZE][SIZE];
	_generateGaussianKernel(m, std);

	for (ull i = distS; i < distE; ++i) {
		x = calcX(i), y = calcY(i);
		p = distPoint(x, y);

		int index = 0, di = i - distS;

		// 处理边界
		if (x < D || x >= distW - D || y < D || y >= distH - D)
			continue;

		fr = fg = fb = 0;
		// 计算核 SIZE*SIZE
		for (int dx = -D; dx <= D; ++dx)
			for (int dy = -D; dy <= D; ++dy) {
				p = distPoint(x + dx, y + dy);
				getColor(tr, tg, tb, p);
				mult = m[dx + D][dy + D];
				fr += tr * mult;
				fg += tg * mult;
				fb += tb * mult;
				if (index++ >= N) break;
			}
		fr = max(min(fr, 255), 0);
		fg = max(min(fg, 255), 0);
		fb = max(min(fb, 255), 0);

		*(distPoint(x, y) + 2) = fr;
		*(distPoint(x, y) + 1) = fg;
		*(distPoint(x, y)) = fb;
	}
}

void ImageProcess::gaussianFilter(ImageFragment * frag) {
	gaussianFilter(frag, param.filter.gausStd);
}

void ImageProcess::wienerFilter(ImageFragment * frag) {

	CImage* dist = frag->dist;
	if (dist == NULL || dist->IsNull()) return;

	const int SIZE = WIENER_KERNEL_SIZE, N = SIZE * SIZE;
	const int D = SIZE >> 1, M = N >> 1;

	int distH = dist->GetHeight(), distW = dist->GetWidth();
	int dPit = dist->GetPitch(), dCnt = dist->GetBPP() / 8;

	double start = frag->start, end = frag->end;

	ull distCnt = distW * distH;
	ull distS = distCnt * start, distE = distCnt * end;
	ull length = distE - distS;

	byte* distData = (byte*)dist->GetBits();

	byte *ps[9], *p; // 临时点
	double pt[3];

	int x, y;

	double mult;
	double m[SIZE][SIZE];

	double noise[3];
	double *mean[3], *variance[3];

	// 初始化
	for (int ch = 0; ch < 3; ++ch) {
		mean[ch] = new double[length];
		variance[ch] = new double[length];
	}

	// 计算噪声 & 方差 & 均值
	for (ull i = distS; i < distE; ++i) {
		x = calcX(i), y = calcY(i);

		// 坐标偏移量
		auto offset = i - distS; // index(x, y, distW) 

		// 处理边界
		if (x < D || x >= distW - D || y < D || y >= distH - D)
			continue;

		int index = 0;

		for (int dx = -D; dx <= D; ++dx)
			for (int dy = -D; dy <= D; ++dy)
				ps[index++] = distPoint(x + dx, y + dy);

		// 三通道
		for (int ch = 0; ch < 3; ++ch) {
			mean[ch][offset] = 0.0;
			variance[ch][offset] = 0.0;
			for (int i = 0; i < N; ++i)
				mean[ch][offset] += ps[i][ch];
			mean[ch][offset] /= N;
			for (int i = 0; i < N; ++i)
				variance[ch][offset] += pow(ps[i][ch] - mean[ch][offset], 2.0);
			variance[ch][offset] /= N;
			noise[ch] += variance[ch][offset];
		}
	}
	for (int ch = 0; ch < 3; ++ch)
		noise[ch] /= length;

	// loop #2: do Wiener filter
	for (ull i = distS; i < distE; ++i) {
		x = calcX(i), y = calcY(i);

		// 坐标偏移量
		auto offset = i - distS; // index(x, y, distW) 

		// 处理边界
		if (x < D || x >= distW - D || y < D || y >= distH - D)
			continue;

		p = distPoint(x, y);

		// 三通道
		for (int ch = 0; ch < 3; ++ch) {
			pt[ch] = p[ch] - mean[ch][offset];
			double t = variance[ch][offset] - noise[ch];
			if (t < 0.0) t = 0.0;
			variance[ch][offset] = fmax(variance[ch][offset], noise[ch]);
			pt[ch] = pt[ch] / variance[ch][offset] * t + mean[ch][offset];

			*(distPoint(x, y) + ch) = pt[ch];
		}
	}
	for (int ch = 0; ch < 3; ++ch) {
		delete[] mean[ch];
		delete[] variance[ch];
	}

}

void ImageProcess::bilateralFiter(ImageFragment * frag, double std) {

	CImage* dist = frag->dist;
	if (dist == NULL || dist->IsNull()) return;

	const int SIZE = BILATERAL_KERNEL_SIZE, N = SIZE * SIZE;
	const int D = SIZE >> 1, M = N >> 1;

	int distH = dist->GetHeight(), distW = dist->GetWidth();
	int dPit = dist->GetPitch(), dCnt = dist->GetBPP() / 8;

	double start = frag->start, end = frag->end;

	ull distCnt = distW * distH;
	ull distS = distCnt * start, distE = distCnt * end;
	ull length = distE - distS;

	byte* distData = (byte*)dist->GetBits();

	byte *p, *cp; // 点, 中心点

	double* mr = new double[N];
	double* mg = new double[N];
	double* mb = new double[N];

	int x, y;
	double fr, fg, fb;
	double tr, tg, tb;
	int dr, dg, db;

	double wr, wg, wb; // 权值
	double sr, sg, sb; // sum

	double mult;
	double m[SIZE][SIZE], c[256];
	_generateGaussianKernelForBin(m, std);
	_generateColorKernel(c, std);

	for (ull i = distS; i < distE; ++i) {
		x = calcX(i), y = calcY(i);
		cp = distPoint(x, y);
		getColor(fr, fg, fb, cp);

		int index = 0, di = i - distS;

		// 处理边界
		if (x < D || x >= distW - D || y < D || y >= distH - D)
			continue;

		sr = sg = sb = 0;
		// 计算核 SIZE*SIZE
		for (int dx = -D; dx <= D; ++dx)
			for (int dy = -D; dy <= D; ++dy) {
				p = distPoint(x + dx, y + dy);
				getColor(tr, tg, tb, p);

				dr = abs(tr - fr); dg = abs(tg - fg); db = abs(tb - fb);

				mult = m[dx + D][dy + D];
				mr[index] = c[dr] * mult;
				mg[index] = c[dg] * mult;
				mb[index] = c[db] * mult;

				sr += mr[index]; sg += mg[index]; sb += mb[index];

				if (index++ >= N) break;
			}

		// 归一化
		for (int j = 0; j < N; ++j) {
			mr[j] /= sr; mg[j] /= sg; mb[j] /= sb;
		}

		index = 0;

		fr = fg = fb = 0;
		for (int dx = -D; dx <= D; ++dx)
			for (int dy = -D; dy <= D; ++dy) {
				p = distPoint(x + dx, y + dy);
				getColor(tr, tg, tb, p);

				fr += tr * mr[index]; // 归一化
				fg += tg * mg[index];
				fb += tb * mb[index];

				if (index++ >= N) break;
			}

		fr = max(min(fr, 255), 0);
		fg = max(min(fg, 255), 0);
		fb = max(min(fb, 255), 0);

		*(distPoint(x, y) + 2) = fr;
		*(distPoint(x, y) + 1) = fg;
		*(distPoint(x, y)) = fb;
	}
}

void ImageProcess::bilateralFiter(ImageFragment * frag) {
	bilateralFiter(frag, param.filter.gausStd);
}

void ImageProcess::saltNoise(ImageFragment * frag, double factor) {
	CImage* dist = frag->dist;
	if (dist == NULL || dist->IsNull()) return;

	int distH = dist->GetHeight(), distW = dist->GetWidth();
	int dPit = dist->GetPitch(), dCnt = dist->GetBPP() / 8;

	double start = frag->start, end = frag->end;

	ull distCnt = distW * distH;
	ull distS = distCnt * start, distE = distCnt * end;

	byte* distData = (byte*)dist->GetBits();

	int x, y;

	for (ull i = distS; i < distE; ++i) {
		x = calcX(i), y = calcY(i);

		if ((rand() / (double)RAND_MAX) <= factor) {
			byte val = (rand() & 0x1) ? 0 : 255;

			*(distPoint(x, y) + 2) = val;
			*(distPoint(x, y) + 1) = val;
			*(distPoint(x, y)) = val;
		}
	}
}

void ImageProcess::saltNoise(ImageFragment * frag) {
	saltNoise(frag, param.noise.factor);
}

double ImageProcess::_generateBoxMuller(double avg, double std) {
	double u1, u2;
	static double z0, z1;
	static bool generated = false;
	generated = !generated;
	if (!generated) return z1 * std + avg;
	do {
		u1 = (double)rand() / RAND_MAX;
		u2 = (double)rand() / RAND_MAX;
	} while (u1 <= DBL_MIN);
	z0 = sqrt(-2.0 * log(u1)) * cos(PI2 * u2);
	z1 = sqrt(-2.0 * log(u1)) * sin(PI2 * u2);
	return z0 * std + avg;
}

void ImageProcess::_generateGaussianKernel(double t[3][3], double std) {
	// std => sigma

	const int SIZE = GAUS_KERNEL_SIZE, N = SIZE * SIZE;
	const int D = SIZE >> 1, M = N >> 1;
	const double f = 1 / (2.0 * PI * std); // 1 / (2.0 * PI * std * std) 
// [[0, 1, 2], [0, 1, 2], [0, 1, 2]], center is (1, 1)
	double total = 0;
	for (int i = 0; i < SIZE; ++i) {
		double xsq = pow(i - D, 2.0);
		for (int j = 0; j < SIZE; ++j) {
			double ysq = pow(j - D, 2.0);
			double e = exp(-(xsq + ysq) / (2.0 * std * std));
			t[i][j] = f * e;
			total += t[i][j];
		}
	}
	for (int i = 0; i < SIZE; ++i)
		for (int j = 0; j < SIZE; ++j)
			t[i][j] /= total;
}

// 用于双边滤波器的高斯核
void ImageProcess::_generateGaussianKernelForBin(double t[3][3], double std) {

	const int SIZE = GAUS_KERNEL_SIZE, N = SIZE * SIZE;
	const int D = SIZE >> 1, M = N >> 1;
// [[0, 1, 2], [0, 1, 2], [0, 1, 2]], center is (1, 1)
	double total = 0;
	for (int i = 0; i < SIZE; ++i) {
		double xsq = pow(i - D, 2.0);
		for (int j = 0; j < SIZE; ++j) {
			double ysq = pow(j - D, 2.0);
			double e = exp(-(xsq + ysq) / (2.0 * std * std));
			t[i][j] = e;
			total += t[i][j];
		}
	}
}

void ImageProcess::_generateColorKernel(double c[256], double std) {
	for (int i = 0; i < 256; ++i)
		c[i] = exp(-(i*i) / (2 * std * std));
}

double ImageProcess::_biCubicPoly(double x) {
	double absX = abs(x);
	double a = -0.5;
	if (absX <= 1.0)
		return (a + 2)*pow(absX, 3) - (a + 3)*pow(absX, 2) + 1;
	else if (absX < 2.0)
		return a * pow(absX, 3) - 5 * a*pow(absX, 2) + 8 * a*absX - 4 * a;
	else
		return 0.0;
	/*
	if (bcCache.find(absX) == bcCache.end()) {
		double a = -0.5;
		if (absX <= 1.0)
			bcCache[absX] = (a + 2)*pow(absX, 3) - (a + 3)*pow(absX, 2) + 1;
		else if (absX < 2.0)
			bcCache[absX] = a * pow(absX, 3) - 5 * a*pow(absX, 2) + 8 * a*absX - 4 * a;
		else
			bcCache[absX] = 0.0;
		LOG("bcCache[" << absX << "] = " << bcCache[absX]);
	}
	return bcCache[absX];
	*/
}

void ImageProcess::_biCubicTranslate(ull distS, ull distE, int distW, int distH,
	int distW2, int distH2, float fCosa, float fSina, double mx, double my,
	byte *srcData, int pit, int cnt, int srcW, int srcH, byte * distData, int dPit, int dCnt) {

	double dx, dy, tx, ty;

	// W(x - xi) W(y - yj)
	double wx, wy;

	int x, y; // 实际坐标
	int minX, maxX, minY, maxY;

	int fr, fb, fg; // 最终颜色
	int tr, tb, tg; // 中间运算

	byte *pi; // 临时点

	for (ull i = distS; i < distE; ++i) {
		int distX = calcX(i), distY = calcY(i);

		dx = distX - distW / 2, dy = -distY + distH / 2;

		// tx = xcosa + ysina, ty = ycosa - xsina
		// tx = tx + distW2/2, ty = -ty + distH2/2 
		tx = dx * fCosa + dy * fSina + distW2 / 2;
		ty = dx * fSina - dy * fCosa + distH2 / 2;

		x = tx / mx, y = ty / my;
		/*
		LOG("dist: " << distX << ", " << distY << " => d:" << dx << ", " << dy << 
			" => t:" << tx << ", " << ty << " => src: " << x << ", " << y);
		*/
		minX = x - 1; maxX = x + 2;
		minY = y - 1; maxY = y + 2;

		fr = fg = fb = 0;

		if ((minX >= 0) && (maxX < srcW) && (minY >= 0) && (maxY < srcH))
			// 在范围内
			for (int xi = minX; xi <= maxX; ++xi)
				for (int yj = minY; yj <= maxY; ++yj) {
					pi = srcPoint(xi, yj);
					wx = _biCubicPoly(x - xi);
					wy = _biCubicPoly(y - yj);
					getColor(tr, tg, tb, pi);
					fr += tr * wx * wy;
					fg += tg * wx * wy;
					fb += tb * wx * wy;
				} else fr = fg = fb = 255;

		*(distPoint(distX, distY) + 2) = fr;
		*(distPoint(distX, distY) + 1) = fg;
		*(distPoint(distX, distY)) = fb;
	}
}

void ImageProcess::_normalTranslate(ull distS, ull distE, int distW, int distH,
	int distW2, int distH2, float fCosa, float fSina, double mx, double my,
	byte *srcData, int pit, int cnt, int srcW, int srcH, byte * distData, int dPit, int dCnt) {

	double dx, dy, tx, ty;

	// W(x - xi) W(y - yj)
	double wx, wy;

	int x, y; // 实际坐标

	// 四个临近点坐标
	int x1, x2, y1, y2;
	byte *p1, *p2, *p3, *p4;

	int fr, fb, fg; // 最终颜色
	int tr1, tb1, tg1; // 中间运算
	int tr2, tb2, tg2; // 中间运算
	double  epsilon = 0.001;

	for (ull i = distS; i < distE; ++i) {
		int distX = calcX(i), distY = calcY(i);

		dx = distX - distW / 2, dy = -distY + distH / 2;

		// tx = xcosa + ysina, ty = ycosa - xsina
		// tx = tx + distW2/2, ty = -ty + distH2/2 
		tx = dx * fCosa + dy * fSina + distW2 / 2;
		ty = dx * fSina - dy * fCosa + distH2 / 2;

		x = tx / mx, y = ty / my;

		fr = fg = fb = 0;

		//计算四个最临近象素的坐标，+1向右下方移动  
		x1 = (int)x; x2 = x1 + 1;
		y1 = (int)y; y2 = y1 + 1;

		if ((x1 >= 0) && (x2 < srcW) && (y1 >= 0) && (y2 < srcH)) {
			p1 = srcPoint(x1, y1); p2 = srcPoint(x2, y1);
			p3 = srcPoint(x1, y2); p4 = srcPoint(x2, y2);

			if (fabs(x - srcW + 1) <= epsilon) {
				// 要计算的点在图像右边缘上  
				if (fabs(y - srcH + 1) <= epsilon)
					// 要计算的点正好是图像最右下角那一个象素，直接返回该点象素值  
					getColor(fr, fg, fb, p1);
				else {
					// 在图像右边缘上且不是最后一点，直接一次插值即可  
					getColor(fr, fg, fb, p1);
					getColor(tr1, tg1, tb1, p3);
					fr = (int)(fr + (y - y1) * (tr1 - fr));
					fg = (int)(fg + (y - y1) * (tg1 - fg));
					fb = (int)(fb + (y - y1) * (tb1 - fb));
				}
			} else if (fabs(y - srcH + 1) <= epsilon) {
				// 要计算的点在图像下边缘上且不是最后一点，直接一次插值即可  
				getColor(fr, fg, fb, p1);
				getColor(tr1, tg1, tb1, p2);
				fr = (int)(fr + (x - x1) * (tr1 - fr));
				fg = (int)(fg + (x - x1) * (tg1 - fg));
				fb = (int)(fb + (x - x1) * (tb1 - fb));
			} else {
				getColor(fr, fg, fb, p1);
				getColor(tr1, tg1, tb1, p2);
				fr = (int)(fr + (x - x1) * (tr1 - fr));
				fg = (int)(fg + (x - x1) * (tg1 - fg));
				fb = (int)(fb + (x - x1) * (tb1 - fb));
				getColor(tr1, tg1, tb1, p3);
				getColor(tr2, tg2, tb2, p4);
				tr1 = (int)(tr1 + (x - x1) * (tr2 - tr1));
				tg1 = (int)(tg1 + (x - x1) * (tg2 - tg1));
				tb1 = (int)(tb1 + (x - x1) * (tb2 - tb1));

				fr = (int)(fr + (y - y1) * (tr1 - fr));
				fg = (int)(fg + (y - y1) * (tg1 - fg));
				fb = (int)(fb + (y - y1) * (tb1 - fb));
			}
		} else fr = fg = fb = 255;

		*(distPoint(distX, distY) + 2) = fr;
		*(distPoint(distX, distY) + 1) = fg;
		*(distPoint(distX, distY)) = fb;
	}

}

/*
void ImageProcess::rotateImage(ImageFragment *frag, double angle) {
	CImage *src = frag->src, *dist = frag->dist;
	if (src == NULL || dist == NULL) return;
	if (src->IsNull() || dist->IsNull()) return;

	double start = frag->start, end = frag->end;

	int srcH = src->GetHeight(), srcW = src->GetWidth();
	int distH = dist->GetHeight(), distW = dist->GetWidth();
	int pit = src->GetPitch(), cnt = src->GetBPP() / 8;
	int dPit = dist->GetPitch(), dCnt = dist->GetBPP() / 8;

	ull srcCnt = srcW * srcH, distCnt = distW * distH;
	ull srcS = srcCnt * start, srcE = srcCnt * end;
	ull distS = distCnt * start, distE = distCnt * end;

	byte* srcData = (byte*)src->GetBits();
	byte* distData = (byte*)dist->GetBits();
	byte* tmpData = frag->tmpData;

	// 临时点
	int x, y; byte *p, *pi;
	// 结果
	int fr, fg, fb;

	float fCosa = cos(angle * PI / 180), fSina = sin(angle * PI / 180);

	float const1 = (float)(-0.5*(srcW - 1)*fCosa - 0.5*(srcH - 1)*fSina + 0.5*(srcW - 1));
	float const2 = (float)(0.5*(srcW - 1)*fSina - 0.5*(srcH - 1)*fCosa + 0.5*(srcH - 1));

	for (ull i = distS; i < distE; ++i) {
		int distX = calcX(i), distY = calcY(i);

		x = (float)distX*fCosa + (float)distY*fSina + const1 + 0.5;
		y = -(float)distX*fSina + (float)distY*fCosa + const2 + 0.5;

		ull id = index(x, y, srcW) / srcCnt;

		if (y >= 0 && y + 2 <= srcH && x >= 0 && x + 2 <= srcW) {
			p = srcPoint(x, y); getColor(fr, fg, fb, p);
		} else
			fr = 255, fg = 255, fb = 255;

		*(tmpPoint(distX, distY) + 2) = fr;
		*(tmpPoint(distX, distY) + 1) = fg;
		*(tmpPoint(distX, distY)) = fb;
	}
}
*/
/*
void ImageProcess::rotateImage(ImageFragment *frag) {
	rotateImage(frag, param.angle);
}
*/
/*
void ImageProcess::scaleImage(ImageFragment *frag, double xScale, double yScale, int mode) {
	CImage *src = frag->src, *dist = frag->dist;
	if (src == NULL || dist == NULL) return;
	switch (mode) {
	case 0:
		_biCubicScale(frag, xScale, yScale); break;
	case 1:
		_normalScale(frag, xScale, yScale); break;
	//case 2:
	//	_2DArrayScale(frag, xScale, yScale); break;
	}
}
*/
/*
void ImageProcess::scaleImage(ImageFragment *frag) {
	scaleImage(frag, param.xScale, param.yScale, param.algo);
}
*/

/*
void ImageProcess::_biCubicScale(ImageFragment * frag, double xScale, double yScale) {
	CImage *src = frag->src, *dist = frag->dist;

	int srcH = src->GetHeight(), srcW = src->GetWidth();
	int distH = dist->GetHeight(), distW = dist->GetWidth();

	int distH2 = srcH * yScale, distW2 = srcW * xScale; // 实际图片所占大小
	int dx = (distW - distW2) / 2, dy = (distH - distH2) / 2; // 偏移

	int pit = src->GetPitch(), cnt = src->GetBPP() / 8;

	double start = frag->start, end = frag->end;

	ull srcCnt = srcW * srcH, distCnt = distW * distH;
	ull srcS = srcCnt * start, srcE = srcCnt * end;
	ull distS = distCnt * start, distE = distCnt * end;

	LOG("from " << start << " to " << end);
	LOG("srcS = " << srcS << ", srcE = " << srcE);
	LOG("distS = " << distS << ", distE = " << distE);
	LOG("srcCnt = " << srcCnt << ", distCnt = " << distCnt);
	LOG("distCnt = " << distW << ", distH = " << distH);

	if (!dist->IsNull() && (dist->GetWidth() != distW || dist->GetHeight() != distH))
		dist->Destroy();
	if (dist->IsNull()) dist->Create(distW, distH, cnt * 8, 0);

	if (distW == 0 || distH == 0) return;

	byte* srcData = (byte*)src->GetBits();
	byte* distData = (byte*)dist->GetBits();

	int dPit = dist->GetPitch(), dCnt = dist->GetBPP() / 8;

	// W(x - xi) W(y - yj)
	double wx, wy;
	byte *p, *pi;

	int fr, fb, fg; // 最终颜色
	int tr, tb, tg; // 中间运算

	int minX, maxX, minY, maxY;

	double mx = ((double)distW2 - 1) / ((double)srcW - 1);
	double my = ((double)distH2 - 1) / ((double)srcH - 1);

	for (ull i = distS; i < distE; ++i) {
		int distX = calcX(i), distY = calcY(i);

		int x = (double)(distX-dx) / mx;
		int y = (double)(distY-dy) / my;
		p = srcPoint(x, y);

		minX = x - 1; maxX = x + 2;
		minY = y - 1; maxY = y + 2;

		fr = fg = fb = 0;

		if ((minX >= 0) && (maxX < srcW) && (minY >= 0) && (maxY < srcH)) 
			for (int xi = minX; xi <= maxX; ++xi)
				for (int yj = minY; yj <= maxY; ++yj) {
					pi = srcPoint(xi, yj);
					wx = _biCubicPoly(x - xi);
					wy = _biCubicPoly(y - yj);
					getColor(tr, tg, tb, pi);
					fr += tr * wx * wy;
					fg += tg * wx * wy;
					fb += tb * wx * wy;
				}
		else fr = fg = fb = 255;

		*(distPoint(distX, distY) + 2) = fr;
		*(distPoint(distX, distY) + 1) = fg;
		*(distPoint(distX, distY)) = fb;
	}
}
*/
/*
void ImageProcess::_normalScale(ImageFragment *frag, double xScale, double yScale) {
	CImage *src = frag->src, *dist = frag->dist;

	int srcH = src->GetHeight(), srcW = src->GetWidth();	
	int distH = dist->GetHeight(), distW = dist->GetWidth();

	int distH2 = srcH * yScale, distW2 = srcW * xScale; // 实际图片所占大小
	int dx = (distW - distW2) / 2, dy = (distH - distH2) / 2; // 偏移

	int pit = src->GetPitch(), cnt = src->GetBPP() / 8;

	double start = frag->start, end = frag->end;

	ull srcCnt = srcW * srcH, distCnt = distW * distH;
	ull srcS = srcCnt * start, srcE = srcCnt * end;
	ull distS = distCnt * start, distE = distCnt * end;

	if (!dist->IsNull() && (dist->GetWidth() != distW || dist->GetHeight() != distH))
		dist->Destroy();
	if (dist->IsNull()) dist->Create(distW, distH, cnt * 8, 0);

	byte* srcData = (byte*)src->GetBits();
	byte* distData = (byte*)dist->GetBits();

	int dPit = dist->GetPitch(), dCnt = dist->GetBPP() / 8;

	// 四个临近点坐标
	int x1, x2, y1, y2;
	byte *p1, *p2, *p3, *p4;

	int fr, fb, fg; // 最终颜色
	int tr1, tb1, tg1; // 中间运算
	int tr2, tb2, tg2; // 中间运算
	double  epsilon = 0.001;

	double mx = ((double)distW2 - 1) / ((double)srcW - 1);
	double my = ((double)distH2 - 1) / ((double)srcH - 1);

	for (ull i = distS; i < distE; ++i) {
		int distX = calcX(i), distY = calcY(i);

		int x = (double)(distX - dx) / mx;
		int y = (double)(distY - dy) / my;

		//计算四个最临近象素的坐标，+1向右下方移动  
		x1 = (int)x; x2 = x1 + 1;
		y1 = (int)y; y2 = y1 + 1;

		if ((x1 >= 0) && (x2 < srcW) && (y1 >= 0) && (y2 < srcH)) {
			p1 = srcPoint(x1, y1); p2 = srcPoint(x2, y1);
			p3 = srcPoint(x1, y2); p4 = srcPoint(x2, y2);

			if (fabs(x - srcW + 1) <= epsilon) {
				// 要计算的点在图像右边缘上  
				if (fabs(y - srcH + 1) <= epsilon)
					// 要计算的点正好是图像最右下角那一个象素，直接返回该点象素值  
					getColor(fr, fg, fb, p1);
				else {
					// 在图像右边缘上且不是最后一点，直接一次插值即可  
					getColor(fr, fg, fb, p1);
					getColor(tr1, tg1, tb1, p3);
					fr = (int)(fr + (y - y1) * (tr1 - fr));
					fg = (int)(fg + (y - y1) * (tg1 - fg));
					fb = (int)(fb + (y - y1) * (tb1 - fb));
				}
			} else if (fabs(y - srcH + 1) <= epsilon) {
				// 要计算的点在图像下边缘上且不是最后一点，直接一次插值即可  
				getColor(fr, fg, fb, p1);
				getColor(tr1, tg1, tb1, p2);
				fr = (int)(fr + (x - x1) * (tr1 - fr));
				fg = (int)(fg + (x - x1) * (tg1 - fg));
				fb = (int)(fb + (x - x1) * (tb1 - fb));
			} else {
				getColor(fr, fg, fb, p1);
				getColor(tr1, tg1, tb1, p2);
				fr = (int)(fr + (x - x1) * (tr1 - fr));
				fg = (int)(fg + (x - x1) * (tg1 - fg));
				fb = (int)(fb + (x - x1) * (tb1 - fb));
				getColor(tr1, tg1, tb1, p3);
				getColor(tr2, tg2, tb2, p4);
				tr1 = (int)(tr1 + (x - x1) * (tr2 - tr1));
				tg1 = (int)(tg1 + (x - x1) * (tg2 - tg1));
				tb1 = (int)(tb1 + (x - x1) * (tb2 - tb1));

				fr = (int)(fr + (y - y1) * (tr1 - fr));
				fg = (int)(fg + (y - y1) * (tg1 - fg));
				fb = (int)(fb + (y - y1) * (tb1 - fb));
			}
		} else fr = fg = fb = 255;

		*(distPoint(distX, distY) + 2) = fr;
		*(distPoint(distX, distY) + 1) = fg;
		*(distPoint(distX, distY)) = fb;
	}
}
*/
/*
void ImageProcess::_2DArrayScale(ImageFragment *frag, double xScale, double yScale) {
	CImage *src = frag->temp, *dist = frag->dist;

	int srcH = src->GetHeight(), srcW = src->GetWidth();
	int distH = srcH * yScale, distW = srcW * xScale;
	int cnt = src->GetBPP();

	if (frag->index <= 0) {
		dist->Destroy();
		dist->Create(distW, distH, cnt * 8, 0);
	}
	// 四个最临近象素的坐标
	int x1, x2, y1, y2;
	// 四个最临近象素值  
	unsigned char f1, f2, f3, f4;
	// 二个插值中间值  
	unsigned char f12, f34;
	//计算结果  
	int fr, fb, fg;
	double  epsilon = 0.001;
	COLORREF pixel11, pixel12, pixel21, pixel22;
	double mx = ((double)distW - 1) / ((double)srcW - 1);
	double my = ((double)distH - 1) / ((double)srcH - 1);
	for (int i = 0; i < distW; i++) {
		for (int j = 0; j < distH; j++) {
			double x = double((double)i / mx);
			double y = double((double)j / my);
			//计算四个最临近象素的坐标，+1向右下方移动  
			x1 = (int)x; x2 = x1 + 1;
			y1 = (int)y; y2 = y1 + 1;
			if ((x < 0) || (x > srcW - 1) || (y < 0) || (y > srcH - 1)) {
				//要计算的点不在源图范围内，返回-1  
				continue;
			} else {
				if (fabs(x - srcW + 1) <= epsilon) {
					// 要计算的点在图像右边缘上  
					if (fabs(y - srcH + 1) <= epsilon) {
						// 要计算的点正好是图像最右下角那一个象素，直接返回该点象素值  
						pixel11 = src->GetPixel(x1, y1);
						f1 = (unsigned char)GetRValue(pixel11);
						fr = (int)f1;
						f1 = (unsigned char)GetGValue(pixel11);
						fg = (int)f1;
						f1 = (unsigned char)GetBValue(pixel11);
						fb = (int)f1;
					} else {
						// 在图像右边缘上且不是最后一点，直接一次插值即可  
						pixel11 = src->GetPixel(x1, y1);
						pixel12 = src->GetPixel(x1, y2);
						f1 = (unsigned char)GetRValue(pixel11);
						f3 = (unsigned char)GetRValue(pixel12);
						fr = (int)(f1 + (y - y1) * (f3 - f1));
						f1 = (unsigned char)GetGValue(pixel11);
						f3 = (unsigned char)GetGValue(pixel12);
						fg = (int)(f1 + (y - y1) * (f3 - f1));
						f1 = (unsigned char)GetBValue(pixel11);
						f3 = (unsigned char)GetBValue(pixel12);
						fb = (int)(f1 + (y - y1) * (f3 - f1));
					}
				} else if (fabs(y - srcH + 1) <= epsilon) {
					// 要计算的点在图像下边缘上且不是最后一点，直接一次插值即可  
					pixel11 = src->GetPixel(x1, y1);
					pixel21 = src->GetPixel(x2, y1);
					f1 = (unsigned char)GetRValue(pixel11);
					f2 = (unsigned char)GetRValue(pixel21);
					fr = (int)(f1 + (x - x1) * (f2 - f1));
					f1 = (unsigned char)GetGValue(pixel11);
					f2 = (unsigned char)GetGValue(pixel21);
					fg = (int)(f1 + (x - x1) * (f2 - f1));
					f1 = (unsigned char)GetBValue(pixel11);
					f2 = (unsigned char)GetBValue(pixel21);
					fb = (int)(f1 + (x - x1) * (f2 - f1));
				} else {
					pixel11 = src->GetPixel(x1, y1);
					pixel12 = src->GetPixel(x1, y2);
					pixel21 = src->GetPixel(x2, y1);
					pixel22 = src->GetPixel(x2, y2);
					// 计算四个最临近象素值  
					f1 = (unsigned char)GetRValue(pixel11);
					f2 = (unsigned char)GetRValue(pixel21);
					f3 = (unsigned char)GetRValue(pixel12);
					f4 = (unsigned char)GetRValue(pixel22);
					f12 = (unsigned char)(f1 + (x - x1) * (f2 - f1));
					f34 = (unsigned char)(f3 + (x - x1) * (f4 - f3));
					fr = (int)(f12 + (y - y1) * (f34 - f12));
					f1 = (unsigned char)GetGValue(pixel11);
					f2 = (unsigned char)GetGValue(pixel21);
					f3 = (unsigned char)GetGValue(pixel12);
					f4 = (unsigned char)GetGValue(pixel22);
					f12 = (unsigned char)(f1 + (x - x1) * (f2 - f1));
					f34 = (unsigned char)(f3 + (x - x1) * (f4 - f3));
					fg = (int)(f12 + (y - y1) * (f34 - f12));
					f1 = (unsigned char)GetBValue(pixel11);
					f2 = (unsigned char)GetBValue(pixel21);
					f3 = (unsigned char)GetBValue(pixel12);
					f4 = (unsigned char)GetBValue(pixel22);
					f12 = (unsigned char)(f1 + (x - x1) * (f2 - f1));
					f34 = (unsigned char)(f3 + (x - x1) * (f4 - f3));
					fb = (int)(f12 + (y - y1) * (f34 - f12));
				}
			}
			dist->SetPixel(i, j, RGB(fr, fg, fb));
		}
	}
}
*/