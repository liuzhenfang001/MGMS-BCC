/**
Copyright (C) 2023/03/08 by Zhenfang Liu, Jianxiong Ye.
Code must not be used, distributed, without written consent by the authors.
*/
#include "head/naive.hpp"
#include "head/TraceStats.hpp"
#include <vector>
#include <algorithm>
#include <cstring>
#include "head/batchdistance.hpp"
#define COUNTMG
namespace Query4 {

/*template<typename bit_t, uint64_t width>
struct BatchBit {
	static const size_t TYPE_BITS_COUNT = sizeof(bit_t) * 8;
	bit_t data[width];

	BatchBit() : data() {
	}
	void negate() {
		for (unsigned i = 0; i < width; ++i) {
			data[i] = ~data[i];
		}
	}

	void setBit(const size_t ix) {
		auto field = ix / TYPE_BITS_COUNT;
		auto field_bit = ix - (field*TYPE_BITS_COUNT);
		data[field] |= BitBaseOp<bit_t>::getSetMask(field_bit);
	}

	bool notZero() {
		bool not_0 = false;
		for (int i = 0; i < width; i++) {
			if (BitBaseOp<bit_t>::notZero(data[i]))
				not_0 = true;
		}
		return not_0;
	}
};

template<typename bit_type>
struct Queueele {
	//uint32_t distance;
	PersonId perId;
	bit_type change;
	Queueele(PersonId perId, bit_type change) : perId(perId), change(change) { };
	Queueele() :perId(0), change(){ };
};*/




//Persons numDiscovered = runRound(subgraph, seen, toVisit, toVisit.size(), personsRemaining);
// XXX: __attribute__((optimize("align-loops")))  this crashed the compiler
Persons __attribute__((hot)) BFSRunner::runRound(const PersonSubgraph& subgraph, Distance* __restrict__ seen, BFSQueue& toVisit, const Persons numToVisit, const Persons numUnseen) {
   Persons numRemainingToVisit=numToVisit;
   Persons numRemainingUnseen=numUnseen;

   do {
      assert(!toVisit.empty());

      const PersonId person = toVisit.front();
      toVisit.pop_front();

      // Iterate over friends
      auto friendsBounds = subgraph.retrieve(person)->bounds();
      while(friendsBounds.first != friendsBounds.second) {
         if (seen[*friendsBounds.first]) {
            ++friendsBounds.first;
            continue;
         }
         toVisit.push_back_pos() = *friendsBounds.first;

         seen[*friendsBounds.first] = true;
         ++friendsBounds.first;
         numRemainingUnseen--;
      }

      numRemainingToVisit--;
   } while(numRemainingToVisit>0 && numRemainingUnseen>0);

   return numUnseen-numRemainingUnseen;
}

}
