
#pragma once

#include <iostream>

#ifdef _DEBUG
#define LOG(x) cout<<x<<endl;
#else
#define LOG(x) ;
#endif
