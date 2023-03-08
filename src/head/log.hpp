/**
Copyright (C) 2023/03/08 by Zhenfang Liu, Jianxiong Ye.
Code must not be used, distributed, without written consent by the authors.
*/
#pragma once

#include <iostream>

namespace tschrono {
   typedef uint64_t Time;

   Time now();

   struct TimeFrame {
      Time startTime;
      Time endTime;
      Time duration;

      void start();
      void end();
   };
}

#ifdef DEBUG
   #ifdef NDBGPRINT
      #define LOG_PRINT(X) std::cout<<tschrono::now()<<" "<<X<<std::endl
   #else
      #define LOG_PRINT(X) std::cerr<<tschrono::now()<<" "<<X<<std::endl
   #endif
#else
   #define LOG_PRINT(X) std::cout<<tschrono::now()<<" "<<X<<std::endl
#endif

#define FATAL_ERROR(X) std::cerr<<"FATAL: "<<X<<std::endl; throw -1
