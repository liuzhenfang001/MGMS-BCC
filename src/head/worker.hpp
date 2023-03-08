/**
Copyright (C) 2023/03/08 by Zhenfang Liu, Jianxiong Ye.
Code must not be used, distributed, without written consent by the authors.
*/
#pragma once

#include "scheduler.hpp"

struct Workers {

   Scheduler scheduler;
   std::vector<std::thread> threads;

   Workers(uint32_t numWorkers);

   void assist(Scheduler& scheduler);
   void close();
};
