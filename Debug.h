
#pragma once

#include <iostream>

#ifdef _DEBUG
#define LOG(x) cout<<x<<endl;
#else
// 为了输出处理时间，还是需要 cout
#define LOG(x) cout<<x<<endl;
#endif
