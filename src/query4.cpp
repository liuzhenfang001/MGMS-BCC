/**
Copyright (C) 2023/03/08 by Zhenfang Liu, Jianxiong Ye.
Code must not be used, distributed, without written consent by the authors.
*/
#include <iostream>
#include <sstream>
#include <random>
//#include <bitset>
#include <cstring>
#include "head/query4.hpp"
#include <cstdlib>

using namespace std;

std::vector<pair<Query4::PersonId,Query4::PersonId>> generateTasks(const uint64_t maxBfs, const Query4::PersonId graphSize, const size_t batchSize) {
   // Determine number of persons
   Query4::PersonId numBfs = graphSize<maxBfs?graphSize:maxBfs;//the total number of BFS

   // Initialize task size as max of minMorselSize and batchSize;
   uint32_t taskSize=batchSize<Query4::minMorselSize?Query4::minMorselSize:batchSize;//the thread number of NFS

   // Increase task size until max morsel tasks limit is met
   while(numBfs/taskSize>Query4::maxMorselTasks) {//Execution times(Task size) cannot exceed maxMorsel Tasks=256000
      taskSize+=batchSize;
   }

   //LOG_PRINT("[TaskGen] Tasks size "<< taskSize<<" maxBFS= "<<maxBfs);

   // Create tasks
   std::vector<pair<Query4::PersonId,Query4::PersonId>> ranges;
   for (unsigned i = 0; i < numBfs; i+=taskSize) {
      Query4::PersonId start = i;
      Query4::PersonId limit = (start+taskSize)>numBfs?numBfs:(start+taskSize);
      ranges.push_back(make_pair(start,limit));
   }

   return ranges;//if taskSize=10,<0,10><10,20><20,30>...<100,numBfs>
}
// Comperator to create correct ordering of results
namespace Query4 {

size_t getMaxMorselBatchSize() {
   char* userBatchSizeStr;
   userBatchSizeStr = getenv("MAX_BATCH_SIZE");
   if(userBatchSizeStr!=nullptr) {
      return atoi(userBatchSizeStr);
   } else {
      return std::numeric_limits<size_t>::max();
   }
}

vector<double> getCloseness(uint32_t* totalPersons,uint64_t* totalDistances,uint32_t totalgraphsize) {
	// new version
	vector<double> res;
	for (int i = 0; i < SAMPLE_NUMS; i++) {
		if (totalDistances[i] > 0 && totalgraphsize > 0 && totalPersons[i] > 0) {
			res.push_back((static_cast<double>(totalPersons[i])  - 1)*(static_cast<double>(totalPersons[i]) - 1) / (static_cast<double>(totalgraphsize)*static_cast<double>(totalDistances[i])));//static_cast<double>((sizeof(BITYPE) * 8 * BITYPE_WIDTH))
		}
		else {
			res.push_back(0.0);
		}
	}
	return res;
	// old version
	/*if (totalDistances > 0 && totalgraphsize > 0 && totalPersons > 0) {
		return ((static_cast<double>(totalPersons)/(sizeof(BITYPE) * 8 * BITYPE_WIDTH) - 1 )*(static_cast<double>(totalPersons) / (sizeof(BITYPE) * 8 * BITYPE_WIDTH) - 1))/ (static_cast<double>(totalgraphsize)*static_cast<double>(totalDistances )/ (sizeof(BITYPE) * 8 * BITYPE_WIDTH));//static_cast<double>((sizeof(BITYPE) * 8 * BITYPE_WIDTH))
	}
	else {
		return 0.0;
	}*/
}


ResultConcatenator::ResultConcatenator(QueryState* state, const char*& resultOut
   #ifdef STATISTICS
   , BatchStatistics& statistics
   #endif
   )
   : state(state), resultOut(resultOut)
   #ifdef STATISTICS
   , statistics(statistics)
   #endif
{
}

void ResultConcatenator::operator()() {
   #ifdef STATISTICS
   statistics.print();
   #endif

   ostringstream output;
   auto& topEntries=state->topResults.getEntries();
   assert(topEntries.size()<=state->k);
   const uint32_t resNum = min(state->k, (uint32_t)topEntries.size());
   for (uint32_t i=0; i<resNum; i++){
      if(i>0) {
         output<<"|";
      }
      output<<state->subgraph.mapInternalNodeId(topEntries[i].first)<<" "<< topEntries[i].second;//Returns the real ID corresponding to the virtual vertex ID
   }
   const auto& outStr = output.str();
   auto resultBuffer = new char[outStr.size()+1];
   outStr.copy(resultBuffer, outStr.size());
   resultBuffer[outStr.size()]=0;
   resultOut = resultBuffer;
   delete state;
}
}
