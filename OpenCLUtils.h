#pragma once

#include <CL/cl.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

static class OpenCLUtils {
public:
	// 初始化函数
	static void initialize();

	static bool getPlatform();
	static bool getDevices();
	static bool getContext();

	// 初始化步骤函数
	static bool createProgram(string fileName, string kernelName);
	static bool createProgram(string filename);

	static bool createCommonQueue();

	// 实际使用函数
	static cl_mem createBuffer(size_t size, void * pointer);
	static bool readBuffer(cl_mem obj, size_t size, void * res);

	static bool setKernelArg(cl_uint index, size_t size, const void * pointer);

	static bool runKernel(cl_uint dim, const size_t * localSize,
		const size_t * globalSize);

	// 垃圾回收
	static void clean();

private:
	static cl_platform_id platform;
	static cl_device_id* devices;
	static cl_context context;
	static cl_kernel kernel;
	static cl_program program;
	static cl_command_queue queue;

	static vector<cl_mem> mems;

	static cl_int ret;

	static int convertToString(string fileName, std::string& s);

	static bool loadProgram(string fileName);
	static bool buildProgram();
	static bool createKernel(string name);
};

