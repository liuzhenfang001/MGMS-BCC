/**
Copyright (C) 2023/03/08 by Zhenfang Liu, Jianxiong Ye.
Code must not be used, distributed, without written consent by the authors.
*/
#include "head/worker.hpp"
#include "head/log.hpp"

Workers::Workers(uint32_t numWorkers) : scheduler() {
   //LOG_PRINT("[Workers] Allocating worker pool with "<< numWorkers << " workers.");
   for (unsigned i = 0; i < numWorkers; ++i) {
      Executor* executor = new Executor(scheduler,i+1, false);
      threads.emplace_back(&Executor::start, executor);
   }
}

void Workers::assist(Scheduler& tasks) {//用LambdaRunner封装一个Lambda表达式，所产生的的Task放在Workers->scheduler->IoTasks
   for (unsigned i = 0; i < threads.size(); ++i) {
      scheduler.schedule(LambdaRunner::createLambdaTask([&tasks, i] {
         Executor(tasks,i+1, false).run();
      }));
   }
}
void Workers::close() {
   scheduler.setCloseOnEmpty();
   for(auto& thread : threads) {
      thread.join();
   }
}
