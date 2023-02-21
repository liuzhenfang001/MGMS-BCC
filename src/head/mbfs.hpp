//Copyright (C) 2014 by Manuel Then, Moritz Kaufmann, Fernando Chirigati, Tuan-Anh Hoang-Vu, Kien Pham, Alfons Kemper, Huy T. Vo
//
//Code must not be used, distributed, without written consent by the authors
#pragma once

#include "TraceStats.hpp"
#include "statistics.hpp"
#include "base.hpp"
#include "batchdistance.hpp"
#include "bitops.hpp"
#include <array>
#include <cstring>
#define COUNTMS
#define SORTED_NEIGHBOR_PROCESSING
//#define DO_PREFETCH
//#define BI_DIRECTIONAl

namespace Query4 {

	// Batch part
	template<typename bit_t, uint64_t width>
	struct BatchBits1 {
		static const size_t TYPE_BITS_COUNT = sizeof(bit_t) * 8;

		bit_t data[width];

		BatchBits1() : data() {
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
	};

	// General HugeBatchBFS loop
	template<typename bit_t = uint64_t, uint64_t width = 1, bool detectSingle = false, uint32_t alpha = 14, uint32_t beta = 24>
	struct MBFSRunner {
		static const size_t TYPE = 3;
		static const size_t WIDTH = width;//2
		static const size_t TYPE_BITS = sizeof(bit_t) * 8;//bit_t=256
#ifdef DO_PREFETCH
		static const unsigned int PREFETCH = 38;
#endif
		static const size_t BATCH_BITS_COUNT = sizeof(bit_t)*width * 8;
		typedef BatchBits1<bit_t, width> Bitset;//表示一个顶点的bfs情况
		static constexpr uint64_t batchSize() {
			return BATCH_BITS_COUNT;
		}

		static void runBatch(std::vector<BatchBFSdata>& bfsData, const PersonSubgraph& subgraph
#ifdef STATISTICS
			, BatchStatistics& statistics
#endif
		) {
#ifdef COUNTMS
			int num_neigber = 0;
			double num_stat = 0;
#endif // COUNTMS

			

			const auto subgraphSize = subgraph.size();
			// Initialize visit lists
			std::array<Bitset*, 2> visitLists;
			for (int a = 0; a < 2; a++) {
				const auto ret = posix_memalign(reinterpret_cast<void**>(&(visitLists[a])), 64, sizeof(Bitset)*subgraphSize);//数据对齐（分配内存首地址，对齐边界，指定分配字节大小）
				if (unlikely(ret != 0)) {
					throw - 1;
				}
				new(visitLists[a]) Bitset[subgraphSize]();
			}


			const uint32_t numQueries = bfsData.size();//本次并行的源点数，即bfs数
			assert(numQueries > 0 && numQueries <= BATCH_BITS_COUNT);

			// Initialize seen vector
			PersonId minPerson = std::numeric_limits<PersonId>::max();
			Bitset* seen;
			const auto ret = posix_memalign(reinterpret_cast<void**>(&seen), 64, sizeof(Bitset)*subgraphSize);
			if (unlikely(ret != 0)) {
				throw - 1;
			}
			new(seen) Bitset[subgraphSize]();

#ifdef BI_DIRECTIONAl
			bool topDown = true;
			uint64_t unexploredEdges = subgraph.numEdges;
			uint64_t visitNeighbors = 0;
			uint32_t frontierSize = numQueries;
#endif

			// Initialize active queries
			for (size_t pos = 0; pos < numQueries; pos++) {
				assert(seen[bfsData[pos].person].isAllZero());
				seen[bfsData[pos].person].setBit(pos);
				assert(!seen[bfsData[pos].person].isAllZero());
				(visitLists[0])[bfsData[pos].person].setBit(pos);
				minPerson = std::min(minPerson, bfsData[pos].person);

#ifdef BI_DIRECTIONAl
				visitNeighbors += subgraph.retrieve(bfsData[pos].person)->size();//计算源点邻居数目和
#endif
			}

			// Initialize iteration workstate
			Bitset processQuery;//processQuery中的每一位，1表示bfs没有搜索完，0表示已经搜索完
			processQuery.negate();//全置为1

			uint32_t queriesToProcess = numQueries;//代表还剩余的bfs数
			alignas(64) uint32_t numDistDiscovered[BATCH_BITS_COUNT];//512，numDistDiscovered每一个数据代表该bfs在一轮搜索中搜索到的新的顶点数
			memset(numDistDiscovered, 0, BATCH_BITS_COUNT * sizeof(uint32_t));

			BatchDistance<bit_t, width> batchDist(numDistDiscovered);

			size_t curToVisitQueue = 0;
			uint32_t nextDistance = 1;//代表bfs的level

			PersonId startPerson = minPerson;

			// Run iterations
			do {
				size_t startTime = tschrono::now();
				Bitset* const toVisit = visitLists[curToVisitQueue];//visit和visitNext轮替
				Bitset* const nextToVisit = visitLists[1 - curToVisitQueue];

				assert(toVisit != nullptr);
				assert(nextToVisit != nullptr);
				bool isempty = true;
				for (int person = 0; person < subgraphSize; person++) {
					bool zero = true;
#ifdef COUNTMS
					num_stat+=1;
#endif // COUNTMS
					for (int i = 0; i< width; i++) {
						if (BitBaseOp<bit_t>::notZero(toVisit[person].data[i])) {
							zero = false;
						}
					}
					if (zero) {
						continue;
					}
					isempty = false;
					const auto& curFriends = *subgraph.retrieve(person);//获得邻居节点
					auto friendsBounds = curFriends.bounds();
					while (friendsBounds.first != friendsBounds.second) {
#ifdef COUNTMS
						num_neigber++;

#endif // COUNTMS
#ifdef COUNTMS
						num_stat+=2;
#endif // COUNTMS
						for (int i = 0; i < width; i++) {

							nextToVisit[*friendsBounds.first].data[i] |= toVisit[person].data[i];
						}
						++friendsBounds.first;
					}
#ifdef COUNTMS
					num_stat +=1;
#endif // COUNTMS
					for (int i = 0; i < width; i++) {
						toVisit[person].data[i] = BitBaseOp<bit_t>::zero();//reset
					}
				}
				for (PersonId curPerson = 0; curPerson < subgraphSize; ++curPerson) {//为什么从0开始？因为nextvisitList记录了邻居，随机分布
					
#ifdef COUNTMS
					
					num_stat += 1;
#endif // COUNTMS
					for (unsigned i = 0; i < width; i++) {
						const bit_t nextVisit = nextToVisit[curPerson].data[i];
						if (BitBaseOp<bit_t>::notZero(nextVisit)) {

#ifdef COUNTMS

							num_stat += 4;
#endif // COUNTMS
							const bit_t newVisits = BitBaseOp<bit_t>::andNot(nextVisit, seen[curPerson].data[i]);
							nextToVisit[curPerson].data[i] = newVisits;
							for (int j = 0; j < TYPE_BITS; j++) {
#ifdef COUNTMS

								num_stat += 1;
#endif // COUNTMS
								if (BitBaseOp<bit_t>::notZero(newVisits & BitBaseOp<bit_t>::getSetMask(j))) {
									bfsData[i*TYPE_BITS+j].totalDistances[0] += nextDistance;
									bfsData[i*TYPE_BITS + j].totalReachable[0] += 1;
								}
							}
							
							if (BitBaseOp<bit_t>::notZero(newVisits)) {

#ifdef COUNTMS

								num_stat += 1;
#endif // COUNTMS
								seen[curPerson].data[i] |= newVisits;
							}
						}
					}

					
				}
//#ifdef COUNTMS
	//			num_stat+= subgraphSize;
//#endif // COUNTMS

				/*for (int a = 0; a < BATCH_BITS_COUNT; a++) {
					batchDist[a].finalize();
					bfsData[a].totalReachable[0] += batchDist[a].numDistDiscovered[0];
					bfsData[a].totalDistances[0] += batchDist[a].numDistDiscovered[0] * nextDistance;
					memset(numDistDiscovered, 0, BATCH_BITS_COUNT * SAMPLE_NUMS * sizeof(uint32_t));
				}*/
				// Update stats for all processed queries and check if the query is finished
#ifdef DEBUG
				uint64_t newReached = 0;
				uint64_t numNotZero = 0;
#endif
				/*for (uint32_t pos = 0; pos < numQueries; pos++) {
					updateProcessQuery(processQuery, pos, numDistDiscovered[pos], bfsData[pos], nextDistance, queriesToProcess);
#ifdef DEBUG
					if (numDistDiscovered[pos] > 0) {
						newReached += numDistDiscovered[pos];
						numNotZero++;
					}
#endif
				}*/

				TraceStats<500>& stats = TraceStats<500>::getStats();
				stats.traceRound(nextDistance);
				stats.addRoundDuration(nextDistance, (tschrono::now() - startTime));

#ifdef COUNTMS
				num_stat /= width;
				stats.addcount(num_neigber, num_stat);
				num_neigber = 0;
				num_stat = 0;
#endif // COUNTMS

				if (isempty) {
					break;
				}
				nextDistance++;

				// As long as not all queries are finished we have to find s.th. every round
				//assert(newReached != 0);

				// Reset iteration state
				//memset(visitLists[curToVisitQueue],0,sizeof(Bitset)*subgraphSize);
				//memset(numDistDiscovered, 0, BATCH_BITS_COUNT * sizeof(uint32_t));//本轮结束，统计下一轮的BFS探索的新定点数

				// Swap queues
				startPerson = 0;
				curToVisitQueue = 1 - curToVisitQueue;
			} while (true);
			//std::cout <<"nextDistance = "<< nextDistance << std::endl;
			free(seen);
			for (int a = 0; a < 2; a++) {
				free(visitLists[a]);
			}
		}


#ifdef SORTED_NEIGHBOR_PROCESSING
		//runBatchRound(subgraph, startPerson, subgraphSize, toVisit, nextToVisit, seen, batchDist, processQuery)
		static std::pair<uint32_t, uint64_t> __attribute__((hot)) runBatchRound(const PersonSubgraph& subgraph, const PersonId startPerson, const PersonId limit, Bitset* visitList, Bitset* nextVisitList, Bitset* __restrict__ seen, BatchDistance<bit_t, width>& batchDist, const Bitset processQuery
#if defined(STATISTICS)
			, BatchStatistics& statistics, uint32_t nextDistance
#elif defined(TRACE)
			, uint32_t nextDistance
#endif
		) {
			// volatile bit_t pref;
#ifdef DO_PREFETCH
			const int p2 = min(PREFETCH, (unsigned int)(limit - startPerson));//PREFETCH=38
			for (int a = 1; a < p2; a++) {
				__builtin_prefetch(visitList + a, 0);//visitList数据预取，不超过38，少了startperson???
				// pref=(visitList + a)->data[0];
			}
#endif

#ifdef TRACE
			size_t numBitsSet[BATCH_BITS_COUNT];
			size_t numFriendsAnalyzed = 0;
			size_t numNewSeen[BATCH_BITS_COUNT];
			memset(numBitsSet, 0, sizeof(size_t)*BATCH_BITS_COUNT);
			memset(numNewSeen, 0, sizeof(size_t)*BATCH_BITS_COUNT);
			std::vector<uint64_t> xx(513);
#endif


			for (PersonId curPerson = startPerson; curPerson < limit; ++curPerson) {
				auto curVisit = visitList[curPerson];

#ifdef DO_PREFETCH
				if (curPerson + PREFETCH < limit) {
					__builtin_prefetch(visitList + curPerson + PREFETCH, 0);
					// pref=(visitList + curPerson + PREFETCH)->data[0];
				}
#endif

#ifdef TRACE
				{
					size_t bitsSet = 0;
					for (int i = 0; i < width; i++) {
						bitsSet += BitBaseOp<bit_t>::popCount(curVisit.data[i]);
					}
					numBitsSet[bitsSet]++;
				}
#endif

				bool zero = true;
				for (int i = 0; i < width; i++) {
					if (BitBaseOp<bit_t>::notZero(curVisit.data[i])) {
						zero = false;
						break;
					}
				}
				if (zero) {
					continue;
				}

				const auto& curFriends = *subgraph.retrieve(curPerson);//获得邻居节点
				auto friendsBounds = curFriends.bounds();
#ifdef DO_PREFETCH
				const int p = min(PREFETCH, (unsigned int)(friendsBounds.second - friendsBounds.first));
				for (int a = 1; a < p; a++) {
					__builtin_prefetch(nextVisitList + *(friendsBounds.first + a), 1);
					// pref=(nextVisitList + *(friendsBounds.first+a))->data[0];
				}
#endif
				//处理
				for (int i = 0; i < width; i++) {
					curVisit.data[i] &= processQuery.data[i];
				}
				while (friendsBounds.first != friendsBounds.second) {
#ifdef DO_PREFETCH
					if (friendsBounds.first + PREFETCH < friendsBounds.second) {
						__builtin_prefetch(nextVisitList + *(friendsBounds.first + PREFETCH), 1);
						// pref=(nextVisitList + *(friendsBounds.first+PREFETCH))->data[0];
					}
#endif

#ifdef TRACE
					{
						numFriendsAnalyzed++;
						size_t newBits = 0;
						for (int i = 0; i < width; i++) {
							newBits += BitBaseOp<bit_t>::popCount(nextVisitList[*friendsBounds.first].data[i] | curVisit.data[i])
								- BitBaseOp<bit_t>::popCount(nextVisitList[*friendsBounds.first].data[i]);
						}
						numNewSeen[newBits]++;
					}
#endif

					for (int i = 0; i < width; i++) {
						nextVisitList[*friendsBounds.first].data[i] |= curVisit.data[i];
					}
					++friendsBounds.first;
				}
				for (int i = 0; i < width; i++) {
					visitList[curPerson].data[i] = BitBaseOp<bit_t>::zero();//reset
				}
			}

#ifdef BI_DIRECTIONAl
			uint32_t frontierSize = 0;
			uint64_t nextVisitNeighbors = 0;
#endif
			for (PersonId curPerson = 0; curPerson < limit; ++curPerson) {//为什么从0开始？因为nextvisitList记录了邻居，随机分布
#ifdef BI_DIRECTIONAl
				bool nextVisitNonzero = false;
#endif
				for (unsigned i = 0; i < width; i++) {
					const bit_t nextVisit = nextVisitList[curPerson].data[i];
					if (BitBaseOp<bit_t>::notZero(nextVisit)) {
						const bit_t newVisits = BitBaseOp<bit_t>::andNot(nextVisit, seen[curPerson].data[i]);
						nextVisitList[curPerson].data[i] = newVisits;
						if (BitBaseOp<bit_t>::notZero(newVisits)) {
							seen[curPerson].data[i] |= newVisits;
							batchDist.updateDiscovered(newVisits, i);

#ifdef BI_DIRECTIONAl
							nextVisitNonzero = true;
#endif
						}
					}
				}
#ifdef BI_DIRECTIONAl
				if (nextVisitNonzero) {
					frontierSize++;//新搜索到的顶点数
					nextVisitNeighbors += subgraph.retrieve(curPerson)->size();//新搜索到的顶点数的邻居数
				}
#endif
			}

/*#ifdef TRACE
			{
				TraceStats<BATCH_BITS_COUNT>& stats = TraceStats<BATCH_BITS_COUNT>::getStats();
				for (int i = 0; i < BATCH_BITS_COUNT; ++i) {
					stats.addRoundVisitBits(nextDistance, i, numBitsSet[i]);
					stats.addRoundFriendBits(nextDistance, i, numNewSeen[i]);
				}
			}
#endif*/

			batchDist.finalize();
#ifdef BI_DIRECTIONAl
			return std::make_pair(frontierSize, nextVisitNeighbors);
#else
			return std::make_pair(0, 0);
#endif
		}

		//frontierInfo = runBatchRoundRev(subgraph, startPerson, subgraphSize, toVisit, nextToVisit, seen, batchDist, processQuery);
		static std::pair<uint32_t, uint64_t> __attribute__((hot)) runBatchRoundRev(const PersonSubgraph& subgraph, const PersonId startPerson, const PersonId limit, Bitset* visitList, Bitset* nextVisitList, Bitset* __restrict__ seen, BatchDistance<bit_t, width>& batchDist, const Bitset/* processQuery*/
#if defined(STATISTICS)
			, BatchStatistics& statistics, uint32_t nextDistance
#elif defined(TRACE)
			, uint32_t nextDistance
#endif
		) {
			uint32_t frontierSize = 0;
			uint64_t nextVisitNeighbors = 0;

			// volatile bit_t pref;
#ifdef DO_PREFETCH
			const int p2 = min(PREFETCH, (unsigned int)(limit - startPerson));//PREFETCH=38
			for (int a = 1; a < p2; a++) {
				__builtin_prefetch(seen + a, 0);//数据预取
			}
#endif

			for (PersonId curPerson = startPerson; curPerson < limit; ++curPerson) {
				auto curSeen = seen[curPerson];

#ifdef DO_PREFETCH
				if (curPerson + PREFETCH < limit) {
					__builtin_prefetch(seen + curPerson + PREFETCH, 0);
				}
#endif

				bool zero = true;
				for (int i = 0; i < width; i++) {
					if (BitBaseOp<bit_t>::notAllOnes(curSeen.data[i])) {
						zero = false;
						break;
					}
				}
				if (zero) {
					continue;
				}

				const auto& curFriends = *subgraph.retrieve(curPerson);
				auto friendsBounds = curFriends.bounds();
#ifdef DO_PREFETCH
				const int p = min(PREFETCH, (unsigned int)(friendsBounds.second - friendsBounds.first));
				for (int a = 1; a < p; a++) {
					__builtin_prefetch(visitList + *(friendsBounds.first + a), 1);
				}
#endif
				Bitset nextVisit;
				while (friendsBounds.first != friendsBounds.second) {
#ifdef DO_PREFETCH
					if (friendsBounds.first + PREFETCH < friendsBounds.second) {
						__builtin_prefetch(visitList + *(friendsBounds.first + PREFETCH), 1);
					}
#endif

					for (int i = 0; i < width; i++) {
						nextVisit.data[i] |= visitList[*friendsBounds.first].data[i];
					}
					++friendsBounds.first;
				}
				for (int i = 0; i < width; i++) {
					nextVisit.data[i] = BitBaseOp<bit_t>::andNot(nextVisit.data[i], curSeen.data[i]);
				}

				nextVisitList[curPerson] = nextVisit;
				bool nextVisitNonzero = false;
				for (int i = 0; i < width; i++) {
					if (BitBaseOp<bit_t>::notZero(nextVisit.data[i])) {
						seen[curPerson].data[i] = curSeen.data[i] | nextVisit.data[i];
						batchDist.updateDiscovered(nextVisit.data[i], i);

						nextVisitNonzero = true;
					}
				}
				if (nextVisitNonzero) {
					frontierSize++;
					nextVisitNeighbors += subgraph.retrieve(curPerson)->size();
				}
			}

			memset(visitList, 0, sizeof(Bitset) * limit);

			// for (PersonId curPerson = 0; curPerson<limit; ++curPerson) {
			//    for(unsigned i=0; i<width; i++) {
			//       const bit_t nextVisit = nextVisitList[curPerson].data[i];
			//       if(BitBaseOp<bit_t>::notZero(nextVisit)) {
			//          const bit_t newVisits = BitBaseOp<bit_t>::andNot(nextVisit, seen[curPerson].data[i]);
			//          nextVisitList[curPerson].data[i] = newVisits;
			//          if(BitBaseOp<bit_t>::notZero(newVisits)) {
			//             seen[curPerson].data[i] |= newVisits;
			//             batchDist.updateDiscovered(newVisits, i);
			//          }
			//       }
			//    }
			// }

			batchDist.finalize();
			return std::make_pair(frontierSize, nextVisitNeighbors);
		}

#else

		enum VisitAction : uint32_t { EMPTY, SINGLE, MULTI };

		struct VisitResult {
			Bitset validVisit;
			uint32_t queryId;
			VisitAction action;
		};

		static void createVisitList(const Bitset& visitList, const Bitset& processQuery, VisitResult& result) {
			if (detectSingle) {
				uint32_t numFields = 0;
				uint32_t nonEmptyField = 0;
				for (unsigned i = 0; i < width; i++) {
					result.validVisit.data[i] = visitList.data[i] & processQuery.data[i];
					if (BitBaseOp<bit_t>::notZero(visitList.data[i])) {
						if (BitBaseOp<bit_t>::notZero(result.validVisit.data[i])) {
							numFields++;
							nonEmptyField = i;
						}
					}
				}

				// Check if only a single bit is set
				if (numFields == 0) {
					result.action = EMPTY;
				}
				else if (numFields == 1) {
					auto bitPos = CtzlOp<bit_t>::ctzl(result.validVisit.data[nonEmptyField]);
					auto masked_field = BitBaseOp<bit_t>::andNot(result.validVisit.data[nonEmptyField], BitBaseOp<bit_t>::getSetMask(bitPos));
					if (BitBaseOp<bit_t>::isZero(masked_field)) {
						result.action = SINGLE;
						result.queryId = nonEmptyField * Bitset::TYPE_BITS_COUNT + bitPos;
					}
					else {
						result.action = MULTI;
					}
				}
				else {
					result.action = MULTI;
				}
			}
			else {
				//Need not detect single
				bool nonZero = false;
				for (unsigned i = 0; i < width; i++) {
					if (BitBaseOp<bit_t>::notZero(visitList.data[i])) {
						nonZero = true;
						break;
					}
					// } else {
					//    result.validVisit.data[i] = BitBaseOp<bit_t>::zero();
					// }
				}
				if (nonZero) {
					for (unsigned i = 0; i < width; i++) {
						result.validVisit.data[i] = visitList.data[i] & processQuery.data[i];
					}
					result.action = MULTI;
				}
				else {
					result.action = EMPTY;
				}
			}
		}

		static void __attribute__((hot)) runBatchRound(const PersonSubgraph& subgraph, const PersonId startPerson, const PersonId limit, Bitset* visitList, Bitset* nextVisitList, Bitset* __restrict__ seen, BatchDistance<bit_t, width>& batchDist, const Bitset processQuery
#if defined(STATISTICS)
			, BatchStatistics& statistics, uint32_t nextDistance
#elif defined(TRACE)
			, uint32_t nextDistance
#endif
		) {
			VisitResult visitResult;

			for (PersonId curPerson = startPerson; curPerson < limit; ++curPerson) {
				createVisitList(visitList[curPerson], processQuery, visitResult);
				// Skip persons with empty visit list
				if (visitResult.action == EMPTY) { continue; }

				const auto& curFriends = *subgraph.retrieve(curPerson);

				//Only single person in this entry
				if (!detectSingle || visitResult.action == MULTI) {
					//More than one person
					auto friendsBounds = curFriends.bounds();
#ifdef DO_PREFETCH
					__builtin_prefetch(seen + *(friendsBounds.first + 1), 0);
					__builtin_prefetch(nextVisitList + *(friendsBounds.first + 1), 1);
					__builtin_prefetch(seen + *(friendsBounds.first + 2), 0);
					__builtin_prefetch(nextVisitList + *(friendsBounds.first + 2), 1);
#endif
					while (friendsBounds.first != friendsBounds.second) {
#ifdef DO_PREFETCH
						if (friendsBounds.first + 3 < friendsBounds.second) {
							__builtin_prefetch(seen + *(friendsBounds.first + 3), 0);
							__builtin_prefetch(nextVisitList + *(friendsBounds.first + 3), 1);
						}
#endif

						for (unsigned i = 0; i < width; i++) {
							bit_t newVisits = BitBaseOp<bit_t>::andNot(visitResult.validVisit.data[i], seen[*friendsBounds.first].data[i]);
							if (BitBaseOp<bit_t>::notZero(newVisits)) {
								seen[*friendsBounds.first].data[i] |= visitResult.validVisit.data[i];
								nextVisitList[*friendsBounds.first].data[i] |= newVisits;

								// Normal case until uint64_t
								batchDist.updateDiscovered(newVisits, i);
							}
						}
						++friendsBounds.first;
					}
				}
				else {
					auto friendsBounds = curFriends.bounds();
					while (friendsBounds.first != friendsBounds.second) {
#ifdef DO_PREFETCH
						if (friendsBounds.first + 3 < friendsBounds.second) {
							__builtin_prefetch(seen + *(friendsBounds.first + 3), 0);
							__builtin_prefetch(nextVisitList + *(friendsBounds.first + 3), 1);
						}
#endif
						auto field = visitResult.queryId / Bitset::TYPE_BITS_COUNT;
						bit_t newVisits = BitBaseOp<bit_t>::andNot(visitResult.validVisit.data[field], seen[*friendsBounds.first].data[field]);
						if (BitBaseOp<bit_t>::notZero(newVisits)) {
							seen[*friendsBounds.first].data[field] |= visitResult.validVisit.data[field];
							nextVisitList[*friendsBounds.first].data[field] |= newVisits;
							batchDist.numDistDiscovered[visitResult.queryId]++;
						}
						++friendsBounds.first;
					}
				}
			}

			batchDist.finalize();
		}

#endif
		//processQuery中的每一位，1表示bfs没有搜索完，0表示已经搜索完
		//updateProcessQuery(processQuery, pos, numDistDiscovered[pos], bfsData[pos], nextDistance, queriesToProcess);
		/*static void updateProcessQuery(Bitset& processQuery, const uint32_t pos, const uint32_t numDiscovered,
			BatchBFSdata& bfsData, const uint32_t distance, uint32_t& queriesToProcess) {
			auto field = pos / Bitset::TYPE_BITS_COUNT;
			auto field_bit = pos - (field*Bitset::TYPE_BITS_COUNT);

			if (BitBaseOp<bit_t>::notZero(processQuery.data[field] & BitBaseOp<bit_t>::getSetMask(field_bit))) {
				bfsData.totalReachable += numDiscovered;//某一个源点总共搜索的顶点，除去源点
				bfsData.totalDistances += numDiscovered * distance;//某一个源点总共搜索的距离

				if ((bfsData.componentSize - 1) == bfsData.totalReachable) {
					processQuery.data[field] = BitBaseOp<bit_t>::andNot(processQuery.data[field], BitBaseOp<bit_t>::getSetMask(field_bit));
					queriesToProcess--;
				}
			}
			else {
				assert(numDiscovered == 0);
			}
		}*/
	};


}