
#pragma once

#include <iostream>

#ifdef _DEBUG
#define LOG(x) cout<<x<<endl;
#else
// Ϊ���������ʱ�䣬������Ҫ cout
#define LOG(x) cout<<x<<endl;
#endif
