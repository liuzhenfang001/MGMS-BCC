//Copyright (C) 2014 by Manuel Then, Moritz Kaufmann, Fernando Chirigati, Tuan-Anh Hoang-Vu, Kien Pham, Alfons Kemper, Huy T. Vo
//
//Code must not be used, distributed, without written consent by the authors
#pragma once

#include "queue.hpp"
#include "graph.hpp"
#include <cstdint>

namespace Query4 {
   typedef uint64_t Distances; // Type for the sum of distances
   typedef uint8_t Distance; // Type for a distance between two points
   typedef uint32_t PersonId; // Type for person ids
   typedef PersonId Persons; // Type for counting persons

   typedef Graph<PersonId> PersonSubgraph;

   inline awfy::FixedSizeQueue<PersonId>& getThreadLocalPersonVisitQueue(size_t queueSize) {
      static __thread awfy::FixedSizeQueue<PersonId>* toVisitPtr=nullptr;//__thread修饰的变量对于各个线程是独立的
      if(toVisitPtr != nullptr) {
         awfy::FixedSizeQueue<PersonId>& q = *toVisitPtr;
         q.reset(queueSize);
         return q;
      } else {
         toVisitPtr = new awfy::FixedSizeQueue<PersonId>(queueSize);
         return *toVisitPtr;
    
      }
   }

   struct BatchBFSdata {
      const PersonId person;
      const Persons componentSize;

	  Distances totalDistances[SAMPLE_NUMS];//探索的距离
	  Persons totalReachable[SAMPLE_NUMS];//探索到的顶点数

      BatchBFSdata(PersonId person, Persons componentSize)
         : person(person), componentSize(componentSize),
           totalDistances(),totalReachable()
      { 
		  memset(totalDistances, 0, sizeof(Distances) * SAMPLE_NUMS);
		  memset(totalReachable, 0, sizeof(Persons) * SAMPLE_NUMS);
	  }

      //BatchBFSdata(const BatchBFSdata&) = delete;
      BatchBFSdata& operator=(const BatchBFSdata&) = delete;

      /*BatchBFSdata(BatchBFSdata&& other)
         : person(other.person), componentSize(other.componentSize),
           totalDistances(other.totalDistances), totalReachable(other.totalReachable)
      { }*/
   };
}