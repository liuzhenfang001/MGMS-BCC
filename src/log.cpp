/**
Copyright (C) 2023/03/08 by Zhenfang Liu, Jianxiong Ye.
Code must not be used, distributed, without written consent by the authors.
*/
#include "head/log.hpp"

#include <chrono>

tschrono::Time tschrono::now() {
   static const std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();
   auto current=std::chrono::high_resolution_clock::now();
   return std::chrono::duration_cast<std::chrono::milliseconds>(current-startTime).count();
}

void tschrono::TimeFrame::start() {
   startTime=tschrono::now();
}

void tschrono::TimeFrame::end() {
   endTime=tschrono::now();
   duration=endTime-startTime;
}
