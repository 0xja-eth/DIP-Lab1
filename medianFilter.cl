#pragma OPENCL EXTENSION cl_amd_printf : enable

#define srcColor(x, y, ch) src[(y) * pitch + (x) * 3 + ch]
__kernel void medianFilterCL(__global uchar *in, int pitch, int width, int height) {
	__global uchar *src = in - pitch * (height - 1);

	int min;
	int window[9];
	int x = get_global_id(0);
	int y = get_global_id(1);

	if (x < width && x > 0 && y < height && y > 0) {

		for(int ch = 0; ch < 3; ++ch) {

			window[0] = (y == 0 || x == 0) ? 0 : srcColor(x-1, y-1, ch);
			window[1] = (y == 0) ? 0 : srcColor(x, y-1, ch);
			window[2] = (y == 0 || x == width - 1) ? 0 : srcColor(x+1, y-1, ch);
			window[3] = (x == 0) ? 0 : srcColor(x-1, y, ch);
			window[4] = srcColor(x, y, ch);
			window[5] = (x == width - 1) ? 0 : srcColor(x+1, y, ch);
			window[6] = (y == height - 1 || x == 0) ? 0 : srcColor(x-1, y+1, ch);
			window[7] = (y == height - 1) ? 0 : srcColor(x, y+1, ch);
			window[8] = (y == height - 1 || x == width - 1) ? 0 : srcColor(x+1, y+1, ch);
			for (unsigned int j = 0; j < 5; j++) {
				min = j;
				for (unsigned int l = j + 1; l < 9; l++)
					if (window[l] < window[min])
						min = l;
				const float temp = window[j];
				window[j] = window[min];
				window[min] = temp;
			}
			srcColor(x, y, ch) = window[4];
		}
	}
}