#pragma OPENCL EXTENSION cl_amd_printf : enable

#define PI (3.14159265358)

#define srcColor(x, y, ch) src[(y) * pit + (x) * 3 + ch]
#define destColor(x, y, ch) dist[(y) * dPit + (x) * 3 + ch]

double _biCubicPoly(double x) {
	double a = -0.5;
	if(x < 0) x = -x;
	if (x < 1.0) return (a + 2.0)*x*x*x - (a + 3.0)*x*x + 1.0;
	else if (x < 2.0) return a * x*x*x - 5.0*a * x*x + 8.0*a * x - 4.0 * a;
	return 0.0;
}

__kernel void transformCL(__global uchar *in, __global uchar *out, int pit, int dPit, int distW, int distH, 
	int srcW, int srcH, double angle, double xScale, double yScale) {
	__global uchar *src = in - pit * (srcH - 1);
	__global uchar *dist = out - dPit * (distH - 1);

	int color[3], tc;

	int distH2 = srcH * yScale, distW2 = srcW * xScale;
	
	double mx = (distW2*1.0 - 1) / (srcW*1.0 - 1);
	double my = (distH2*1.0 - 1) / (srcH*1.0 - 1);
	
	double fCosa = cos(angle * PI / 180), fSina = sin(angle * PI / 180);

	// 循环体	
	int distX = get_global_id(0);
	int distY = get_global_id(1);

	double dx = distX - distW / 2, dy = -distY + distH / 2;

	double tx = dx * fCosa + dy * fSina + distW2 / 2;
	double ty = dx * fSina - dy * fCosa + distH2 / 2;

	int x = tx / mx, y = ty / my;

	int minX = x - 1, maxX = x + 2;
	int minY = y - 1, maxY = y + 2;
	
	double wx, wy;

	for (int ch = 0; ch < 3; ++ch) 
		color[ch] = 0;

	if ((minX >= 0) && (maxX < srcW) && (minY >= 0) && (maxY < srcH))
		// 在范围内
		for (int xi = minX; xi <= maxX; ++xi)
			for (int yj = minY; yj <= maxY; ++yj) {
				wx = _biCubicPoly(x - xi);
				wy = _biCubicPoly(y - yj);

				for (int ch = 0; ch < 3; ++ch) {
					tc = srcColor(x, y, ch);
					color[ch] += tc * wx * wy;
				}
	} else for (int ch = 0; ch < 3; ++ch) 
		color[ch] = 255;

	for (int ch = 0; ch < 3; ++ch) 
		destColor(distX, distY, ch) = color[ch];
}