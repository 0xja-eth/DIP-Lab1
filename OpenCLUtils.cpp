#include "stdafx.h"
#include "OpenCLUtils.h"

#include "Debug.h"

cl_platform_id OpenCLUtils::platform;
cl_device_id* OpenCLUtils::devices;
cl_context OpenCLUtils::context;
cl_kernel OpenCLUtils::kernel;
cl_program OpenCLUtils::program;
cl_command_queue OpenCLUtils::queue;

std::vector<cl_mem> OpenCLUtils::mems;

cl_int OpenCLUtils::ret;

void OpenCLUtils::initialize() {
	ASSERT(getPlatform());
	ASSERT(getDevices());
	ASSERT(getContext());
	LOG("ÕÍ≥…OpenCLº”‘ÿ");
}

bool OpenCLUtils::getPlatform() {
	cl_uint num; platform = 0;
	if (clGetPlatformIDs(0, NULL, &num) != CL_SUCCESS)
		return false;
	if (num <= 0) return false;
	auto platforms = new cl_platform_id[num];
	clGetPlatformIDs(num, platforms, NULL);
	platform = platforms[0];
	delete[] platforms;
	return true;
}

bool OpenCLUtils::getDevices() {
	cl_uint num;
	if (clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 0, NULL, &num) != CL_SUCCESS) return false;
	if (num <= 0) {
		devices = NULL;
		return false;
	}
	devices = new cl_device_id[num];
	clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, num, devices, NULL);
	return true;
}

bool OpenCLUtils::getContext() {
	int status;
	context = clCreateContext(NULL, 1, devices, 
		NULL, NULL, &status);
	return status == CL_SUCCESS;
}

bool OpenCLUtils::createProgram(string fileName, string kernelName) {
	return loadProgram(fileName) && buildProgram() && createKernel(kernelName);
}
bool OpenCLUtils::createProgram(string fileName) {
	return createProgram(fileName, fileName);
}

bool OpenCLUtils::createCommonQueue() {
	queue = clCreateCommandQueue(context, devices[0], 0, &ret);
	return ret == CL_SUCCESS;
}

cl_mem OpenCLUtils::createBuffer(size_t size, void * pointer) {
	cl_mem mem = NULL;
	mem = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &ret);
	if (ret != CL_SUCCESS) return NULL;
	ret = clEnqueueWriteBuffer(queue, mem, CL_TRUE, 0, size, pointer, 0, NULL, NULL);
	if (ret != CL_SUCCESS) return NULL;
	mems.push_back(mem);
	return mem;
}

bool OpenCLUtils::readBuffer(cl_mem obj, size_t size, void * res) {
	ret = clEnqueueReadBuffer(queue, obj, CL_TRUE, 0, size, res, 0, NULL, NULL);
	return ret == CL_SUCCESS;
}

bool OpenCLUtils::setKernelArg(cl_uint index, size_t size, const void * pointer) {
	ret = clSetKernelArg(kernel, index, size, pointer);
	return ret == CL_SUCCESS;
}

bool OpenCLUtils::runKernel(cl_uint dim, const size_t * localSize, const size_t * globalSize) {
	// Launch
	cl_event event;
	ret = clEnqueueNDRangeKernel(
		queue, kernel, dim, NULL, globalSize,
		localSize, 0, NULL, &event
	);
	if (ret != CL_SUCCESS) return false;
	clWaitForEvents(1, &event);
	clReleaseEvent(event);
	return true;
}

void OpenCLUtils::clean() {
	for (auto i : mems) 
		clReleaseMemObject(i);
	mems.clear();
	clReleaseKernel(kernel);
	clReleaseProgram(program);
	clReleaseCommandQueue(queue);
	//clReleaseContext(context);
}

bool OpenCLUtils::loadProgram(string filename) {
	string sourceStr;
	convertToString(filename, sourceStr);
	const char *source = sourceStr.c_str();
	size_t sourceSize[] = { strlen(source) };
	program = clCreateProgramWithSource(context, 1, &source, sourceSize, &ret);
	return ret == CL_SUCCESS;
}

bool OpenCLUtils::buildProgram() {
	ret = clBuildProgram(program, 1, devices, NULL, NULL, NULL);
	return ret == CL_SUCCESS;
}

bool OpenCLUtils::createKernel(string kernelName) {
	kernel = clCreateKernel(program, kernelName.c_str(), &ret);
	return ret == CL_SUCCESS;
}

/** convert the kernel file into a string */
int OpenCLUtils::convertToString(string fileName, string& s) {
	size_t size; char* str;
	std::fstream f(fileName, (std::fstream::in | std::fstream::binary));

	if (f.is_open()) {
		size_t fileSize;
		f.seekg(0, std::fstream::end);
		size = fileSize = (size_t)f.tellg();
		f.seekg(0, std::fstream::beg);
		str = new char[size + 1];
		if (!str) {
			f.close(); return 0;
		}
		f.read(str, fileSize); f.close();
		str[size] = '\0'; s = str;
		delete[] str; return 0;
	}
	cout << "Error: failed to open file\n:" << fileName << endl;
	return -1;
}
