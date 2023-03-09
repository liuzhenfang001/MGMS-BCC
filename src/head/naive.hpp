/**
Copyright (C) 2023/03/08 by Zhenfang Liu, Jianxiong Ye.
Code must not be used, distributed, without written consent by the authors.
*/
#pragma once

#include "base.hpp"
#include "statistics.hpp"
#include <vector>
#include <algorithm>
#include <cstring>
#include "batchdistance.hpp"
#include "TraceStats.hpp"
#include <unistd.h>
#define COUNTMG
namespace Query4 {

	template<typename bit_t, uint64_t width>
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
		Queueele() :perId(0), change() { };
	};


struct BFSRunner {
   static const size_t TYPE=1;
   static const size_t WIDTH=1;
   static const size_t TYPE_BITS=8;

   typedef awfy::FixedSizeQueue<PersonId> BFSQueue;

   static constexpr size_t batchSize() {
      return 1;
   }

   static inline void runBatch(std::vector<BatchBFSdata>& bfsData, const PersonSubgraph& subgraph
      #ifdef STATISTICS
      , BatchStatistics& statistics
      #endif
      ) {
      for(size_t i=0; i<bfsData.size(); i++) {
         run(bfsData[i].person, subgraph, bfsData[i]);
      }
   }

	static void run(const PersonId start, const PersonSubgraph& subgraph, BatchBFSdata& bfsData) {
		const auto subgraphSize = subgraph.size();
		const auto sourcePer = bfsData.person;

		BatchBit<BITYPE, BITYPE_WIDTH>* verstatus;
		{
			const auto ret = posix_memalign(reinterpret_cast<void**>(&verstatus), 128, sizeof(BatchBit<BITYPE, BITYPE_WIDTH>)*subgraphSize);//
			if (unlikely(ret != 0)) {
				std::cout << "unlikely" << std::endl;
				throw - 1;
			}
			new(verstatus) BatchBit<BITYPE, BITYPE_WIDTH>[subgraphSize]();
		}
		BatchBit<BITYPE, BITYPE_WIDTH>* change_verstatus;
		{
			const auto ret = posix_memalign(reinterpret_cast<void**>(&change_verstatus), 128, sizeof(BatchBit<BITYPE, BITYPE_WIDTH>)*subgraphSize);//
			if (unlikely(ret != 0)) {
				std::cout << "unlikely" << std::endl;
				throw - 1;
			}
			new(change_verstatus) BatchBit<BITYPE, BITYPE_WIDTH>[subgraphSize]();
		}

		//size_t maxqueue = 10000000;
		//awfy::FixedSizeQueue<Queueele<BatchBit<BITYPE, BITYPE_WIDTH>>> Change_Q = awfy::FixedSizeQueue<Queueele<BatchBit<BITYPE, BITYPE_WIDTH>>>(maxqueue);

		//awfy::FixedSizeQueue<PersonId> BFS_Queue = awfy::FixedSizeQueue<PersonId>(maxqueue);
		//awfy::FixedSizeQueue<PersonId> next_BFS_Queue = awfy::FixedSizeQueue<PersonId>(maxqueue);
		//BFS_Queue.push_back_pos() = sourcePer;
		verstatus[sourcePer].negate();
		//size_t maxQ = 0;
		//PersonId minPerson = std::numeric_limits<PersonId>::max();

		uint32_t maxvec = subgraphSize;
		//vector<PersonId> VectorPer;
		//vector<PersonId> next_VectorPer;
		//VectorPer.reserve(maxvec);
		//VectorPer.push_back(sourcePer);
		//next_VectorPer.reserve(maxvec);

		alignas(64) uint32_t numDistDiscovered[SAMPLE_NUMS];//512
		memset(numDistDiscovered, 0, SAMPLE_NUMS * sizeof(uint32_t));
		BatchDistance<BITYPE, BITYPE_WIDTH> batchDist(numDistDiscovered);


		alignas(64) uint32_t* vec_flag;
		vec_flag = (uint32_t*)malloc(maxvec * sizeof(uint32_t));

		alignas(64) uint32_t* next_vec_flag;
		next_vec_flag = (uint32_t*)malloc(maxvec * sizeof(uint32_t));

		memset(vec_flag, 0, maxvec * sizeof(uint32_t));
		memset(next_vec_flag, 0, maxvec * sizeof(uint32_t));
		vec_flag[sourcePer] = 1;

		uint32_t cur_vector = 0;
		uint32_t nextDistance = 1;
		uint32_t cur_visit = 0;
		bool nobreak = true;
		//std::cout << "while";

#ifdef COUNTMG
		int num_neigber = 0;
		double num_stat = 0;
#endif // COUNTMG


		while (nobreak) {
			nobreak = false;
			size_t startTime = tschrono::now();
			auto& VEC_flag = cur_vector == 0 ? vec_flag : next_vec_flag;
			auto& next_VEC_flag = cur_vector == 1 ? vec_flag : next_vec_flag;

			//for (int ii = 0; ii < subgraphSize; ii++) {
				//std::bitset<64> bis(verstatus[ii].data[0]);
				//std::bitset<8> bis1(verstatus[ii].data[1]);
				//std::cout << subgraph.mapInternalNodeId(ii) << " " << bis << std::endl;
			//}
			//std::cout << std::endl;
			for (PersonId curPerson = 0; curPerson < subgraphSize; curPerson++) {
				if (VEC_flag[curPerson] == 0)
					continue;
				const auto& curFriends = *subgraph.retrieve(curPerson);
				auto friendsBounds = curFriends.bounds();
				int numfriends = 0;
				int numf = 0;
				while (friendsBounds.first != friendsBounds.second) {
#ifdef COUNTMG
					num_neigber++;
#endif
					BatchBit<BITYPE, BITYPE_WIDTH> D;
					if (VEC_flag[*friendsBounds.first] > 0) {//处理冲突
						bool insert = false;
						for (int s = 0; s < BITYPE_WIDTH; s++) {

#ifdef COUNTMG
							num_stat += 2;
#endif //
							D.data[s] = verstatus[curPerson].data[s] & subgraph.edgesbit[curPerson][numfriends].data[s] & ~verstatus[*friendsBounds.first].data[s];
							if (BitBaseOp<BITYPE>::notZero(D.data[s])) {
#ifdef COUNTMG
								num_stat += 1;
#endif //
								//std::cout << "cruch in distance = "<< nextDistance << std::endl;

								insert = true;
								//uint64_t numofone = BitBaseOp<BITYPE>::popCount(D.data[s]);
								//bfsData.totalDistances += numofone * nextDistance;
								//bfsData.totalReachable += numofone;
							}
						}
						if (insert) {
							//std::cout << subgraph.mapInternalNodeId(curPerson) << "->" << subgraph.mapInternalNodeId(*friendsBounds.first) << std::endl;
							nobreak = true;
							for (int s = 0; s < BITYPE_WIDTH; s++) {
#ifdef COUNTMG
								num_stat += 2;
#endif //
								change_verstatus[*friendsBounds.first].data[s] |= D.data[s];

							}
							numf++;
							//Change_Q.push_back_pos() = Queueele<BatchBit<BITYPE, BITYPE_WIDTH>>(*friendsBounds.first, D);
							//vector<PersonId>::iterator resu_ = find(next_vectorper.begin(), next_vectorper.end(), *friendsBounds.first);
							//next_bfs_queue.push_back_pos() = *friendsBounds.first;
							next_VEC_flag[*friendsBounds.first] = 1;
						}
					}
					else {//No conflicts
						//std::cout << "no cruch in distance = " << nextDistance << std::endl;
						bool insert = false;

						for (int s = 0; s < BITYPE_WIDTH; s++) {
#ifdef COUNTMG
							num_stat += 2;
#endif //

							D.data[s] = verstatus[curPerson].data[s] & subgraph.edgesbit[curPerson][numfriends].data[s] & ~verstatus[*friendsBounds.first].data[s];
							if (BitBaseOp<BITYPE>::notZero(D.data[s])) {
#ifdef COUNTMG
								num_stat += 4;
#endif //
								verstatus[*friendsBounds.first].data[s] |= D.data[s];
								insert = true;
								uint64_t numofone = BitBaseOp<BITYPE>::popCount(D.data[s]);
								batchDist.updateDiscovered(D.data[s], s);
								//bfsData.totalDistances[s] += numofone * nextDistance;
								//bfsData.totalReachable += numofone;
							}
						}
						if (insert) {
							//std::cout << subgraph.mapInternalNodeId(curPerson) << ">>>" << subgraph.mapInternalNodeId(*friendsBounds.first) << std::endl;

							nobreak = true;
							numf++;
							//next_bfs_queue.push_back_pos() = *friendsBounds.first;
							next_VEC_flag[*friendsBounds.first] = 1;
						}
					}
					friendsBounds.first++;
					numfriends++;
				}


				//std::cout << subgraph.mapInternalNodeId(curPerson) << " " << numf << "|";
			}
			for (int per = 0; per < subgraphSize; per++) {
				//int numofone = 0;

				for (int s = 0; s < BITYPE_WIDTH; s++) {

#ifdef COUNTMG
					num_stat += 6;
#endif // COUNTMG
					BITYPE DD = change_verstatus[per].data[s] & ~verstatus[per].data[s];
					verstatus[per].data[s] |= DD;
					batchDist.updateDiscovered(DD, s);
					//numofone += BitBaseOp<BITYPE>::popCount(change_verstatus[per].data[s]);
					change_verstatus[per].data[s] = BitBaseOp<BITYPE>::zero();
				}

				//bfsData.totalDistances += numofone * nextDistance;
				//bfsData.totalReachable += numofone;
			}

			/*for (; !Change_Q.empty(); Change_Q.pop_front()) {
				Queueele<BatchBit<BITYPE, BITYPE_WIDTH>> change_ele = Change_Q.front();
				for (int ptr = 0; ptr < BITYPE_WIDTH; ptr++) {
					//std::cout << "change" << std::endl;
					verstatus[change_ele.perId].data[ptr] |= change_ele.change.data[ptr];
				}
			}
			Change_Q.reset(maxqueue);*/



			for (int l = 0; l < SAMPLE_NUMS; l++) {
				batchDist.finalize();
				bfsData.totalReachable[l] += batchDist.numDistDiscovered[l];
				bfsData.totalDistances[l] += batchDist.numDistDiscovered[l] * nextDistance;
			}
			memset(numDistDiscovered, 0, SAMPLE_NUMS * sizeof(uint32_t));


			TraceStats<500>& stats = TraceStats<500>::getStats();
			stats.traceRound(nextDistance);
			stats.addRoundDuration(nextDistance, (tschrono::now() - startTime));
#ifdef COUNTMG
			num_stat /= 4;
			stats.addcount(num_neigber, num_stat);
			num_neigber = 0;
			num_stat = 0;
#endif // COUNTMG


			memset(VEC_flag, 0, subgraphSize * sizeof(uint32_t));

			cur_vector = 1 - cur_vector;
			nextDistance++;
			/*size_t startTime = tschrono::now();

			auto& bfs_queue = cur_queue == 0 ? BFS_Queue : next_BFS_Queue;
			auto& next_bfs_queue = cur_queue == 1 ? BFS_Queue : next_BFS_Queue;
			//auto& vectorper = cur_vector == 0 ? VectorPer : next_VectorPer;
			//auto& next_vectorper = cur_vector == 1 ? VectorPer : next_VectorPer;
			auto& VEC_flag = cur_vector == 0 ? vec_flag : next_vec_flag;
			auto& next_VEC_flag = cur_vector == 1 ? vec_flag : next_vec_flag;


			for (; !Change_Q.empty(); Change_Q.pop_front()) {
				Queueele<BatchBit<BITYPE, BITYPE_WIDTH>> change_ele = Change_Q.front();
				for (int ptr = 0; ptr < BITYPE_WIDTH; ptr++) {
					//std::cout << "change" << std::endl;
					verstatus[change_ele.perId].data[ptr] |= change_ele.change.data[ptr];
				}
			}
			Change_Q.reset(maxqueue);

			for (int ii = 0; ii < subgraphSize; ii++) {
				std::bitset<8> bis(verstatus[ii].data[0]);
				std::bitset<8> bis1(verstatus[ii].data[1]);
				std::cout << bis << " " << bis1 << std::endl;
			}
			std::cout << std::endl;



			while (!bfs_queue.empty()) {

				PersonId curPerson = bfs_queue.front();
				bfs_queue.pop_front();

				const auto& curFriends = *subgraph.retrieve(curPerson);
				auto friendsBounds = curFriends.bounds();
				int numfriends = 0;
				while (friendsBounds.first != friendsBounds.second ) {
					//PersonId curfriend = *friendsBounds.first;
					BatchBit<BITYPE, BITYPE_WIDTH> D;
					//vector<PersonId>::iterator result = find(vectorper.begin(), vectorper.end(), *friendsBounds.first);

					if (VEC_flag[*friendsBounds.first]>0) {//处理冲突

						bool insert = false;
						for (int s = 0; s < BITYPE_WIDTH; s++) {
							D.data[s] = verstatus[curPerson].data[s] & subgraph.edgesbit[curPerson][numfriends].data[s] & ~verstatus[*friendsBounds.first].data[s];
							if (BitBaseOp<BITYPE>::notZero(D.data[s])) {
								//std::cout << "cruch in distance = "<< nextDistance << std::endl;
								insert = true;
								uint64_t numofone = BitBaseOp<BITYPE>::popCount(D.data[s]);
								bfsData.totalDistances += numofone * nextDistance;
								bfsData.totalReachable += numofone;
							}
						}
						if (insert) {
							Change_Q.push_back_pos() = Queueele<BatchBit<BITYPE, BITYPE_WIDTH>>(*friendsBounds.first, D);
							//vector<PersonId>::iterator resu_ = find(next_vectorper.begin(), next_vectorper.end(), *friendsBounds.first);
							if (next_VEC_flag[*friendsBounds.first]>0) {
								next_bfs_queue.push_back_pos() = *friendsBounds.first;
								next_VEC_flag[*friendsBounds.first]++;
							}
						}
					}
					else {//无冲突
						//std::cout << "no cruch in distance = " << nextDistance << std::endl;
						bool insert = false;
						for (int s = 0; s < BITYPE_WIDTH; s++) {
							D.data[s] = verstatus[curPerson].data[s] & subgraph.edgesbit[curPerson][numfriends].data[s] & ~verstatus[*friendsBounds.first].data[s];
							if (BitBaseOp<BITYPE>::notZero(D.data[s])) {
								verstatus[*friendsBounds.first].data[s] |= D.data[s];
								insert = true;
								uint64_t numofone = BitBaseOp<BITYPE>::popCount(D.data[s]);
								bfsData.totalDistances += numofone * nextDistance;
								bfsData.totalReachable += numofone;
							}
						}
						if (insert) {
							next_bfs_queue.push_back_pos() = *friendsBounds.first;
							next_VEC_flag[*friendsBounds.first]++;
						}
					}
					friendsBounds.first++;
					numfriends++;
				}
			}

			TraceStats<500>& stats = TraceStats<500>::getStats();
			stats.traceRound(nextDistance);
			stats.addRoundDuration(nextDistance, (tschrono::now() - startTime));

			//std::cout << nextDistance << "\t bfsq=" << bfs_queue.size() <<"\t nexbfsq="<< next_bfs_queue.size()<<  "\t chaq="<<Change_Q.size() <<"\t\t vec="<< vectorper.size()<< "\t nextvec=" << next_vectorper.size() << std::endl;

			if (next_bfs_queue.empty() && Change_Q.empty())
				break;

			cur_queue = 1 - cur_queue;
			//vectorper.clear();
			memset(VEC_flag, 0, subgraphSize * sizeof(uint32_t));

			cur_vector = 1 - cur_vector;
			nextDistance++;*/
		}
		//for (int ii = 0; ii < subgraphSize; ii++) {
			//std::bitset<64> bis(verstatus[ii].data[0]);
			//std::bitset<8> bis1(verstatus[ii].data[1]);
			//std::cout << subgraph.mapInternalNodeId(ii) << " " << bis << std::endl;
		//}
		free(verstatus);
		free(change_verstatus);
		/*BFSQueue& toVisit = getThreadLocalPersonVisitQueue(subgraph.size());
		assert(toVisit.empty()); //Data structures are in a sane state

		// Initialize BFS
		Distance* seen = new Distance[subgraph.size()]();
		seen[start] = 1; // Level = distance + 1, Level 0 = not seen
		toVisit.push_back_pos() = start;

		// Run rounds until we can either early exit or have reached all nodes
		Distance distance=0;
		do {
		   size_t startTime = tschrono::now();

		   const Persons personsRemaining=(bfsData.componentSize-1)-bfsData.totalReachable;//
		   Persons numDiscovered = runRound(subgraph, seen, toVisit, toVisit.size(), personsRemaining);
		   assert(distance<std::numeric_limits<Distance>::max());
		   distance++;

		   // Adjust result
		   bfsData.totalReachable+=numDiscovered;
		   bfsData.totalDistances+=numDiscovered*distance;

		   TraceStats<TYPE_BITS*WIDTH>& stats = TraceStats<TYPE_BITS*WIDTH>::getStats();
		   stats.traceRound(distance);
		   stats.addRoundDuration(distance, (tschrono::now()-startTime));

		   // Exit criteria for full BFS
		   if((bfsData.componentSize-1)==bfsData.totalReachable) {
			  break;
		   }
		} while(true);

		delete[] seen;*/

   }

   static Persons runRound(const PersonSubgraph& subgraph, Distance* __restrict__ seen, BFSQueue& toVisit,
      const Persons numToVisit, const Persons numUnseen);
};

}
