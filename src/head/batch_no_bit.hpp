#pragma once

#include "TraceStats.hpp"
#include "statistics.hpp"
#include "base.hpp"
#include "batchdistance.hpp"

#include <array>
#include <cstring>

#define SORTED_NEIGHBOR_PROCESSING
//#define DO_PREFETCH
//#define BI_DIRECTIONAl

namespace Query4 {
	// Batch part
	template<typename bit_t, uint64_t width>
	struct BatchBits {
		static const size_t TYPE_BITS_COUNT = sizeof(bit_t) * 8;

		bit_t data[width];

		BatchBits() : data() {
		}

#ifdef DEBUG
		bool isAllZero() const {
			for (unsigned i = 0; i < width; ++i) {
				if (BitBaseOp<bit_t>::notZero(data[i])) {
					return false;
				}
			}

			return true;
		}
#endif

#ifdef STATISTICS
		size_t count() const {
			size_t sum = 0;
			for (unsigned i = 0; i < width; ++i) {
				sum += BitBaseOp<bit_t>::popCount(data[i]);
			}
			return sum;
		}
#endif

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



	template<typename bit_type>
	struct Quadruple {
		//uint32_t distance;
		PersonId perId;

		uint32_t bfsId;
		bit_type change;
		Quadruple(PersonId perId, bit_type change, uint32_t bfsId) : perId(perId), change(change), bfsId(bfsId) { };
		Quadruple() :perId(), change(), bfsId() { };
	};



	// General HugeBatchBFS loop
	template<typename bit_t = uint64_t, uint64_t width = 1, bool detectSingle = false, uint32_t alpha = 14, uint32_t beta = 24>
	struct HugeBatchBfs {
		static const size_t TYPE = 3;
		static const size_t WIDTH = width;//2
		static const size_t TYPE_BITS = sizeof(bit_t) * 8;//bit_t=256
#ifdef DO_PREFETCH
		static const unsigned int PREFETCH = 38;
#endif
		static const size_t BATCH_BITS_COUNT = sizeof(bit_t)*width * 8;
		typedef BatchBits<bit_t, width> Bitset;//表示一个顶点的bfs情况

		static constexpr uint64_t batchSize() {
			return BATCH_BITS_COUNT;
		}

		static void runBatch(std::vector<BatchBFSdata>& bfsData, const PersonSubgraph& subgraph
#ifdef STATISTICS
			, BatchStatistics& statistics
#endif
		) {
			const uint32_t numQueries = bfsData.size();//本次并行的源点数，即bfs数
			assert(numQueries > 0 && numQueries <= BATCH_BITS_COUNT);
			const auto subgraphSize = subgraph.size();


			BatchBits<BITYPE, BITYPE_WIDTH>** verstatus;
			verstatus = new BatchBits<BITYPE, BITYPE_WIDTH>*[subgraphSize]();
			for (int a = 0; a < subgraphSize; a++) {
				const auto ret = posix_memalign(reinterpret_cast<void**>(&(verstatus[a])), 64, sizeof(BatchBits<BITYPE, BITYPE_WIDTH>)*numQueries);//数据对齐（分配内存首地址，对齐边界，指定分配字节大小）
				if (unlikely(ret != 0)) {
					std::cout << "unlikely" << std::endl;
					throw - 1;
				}
				new(verstatus[a]) BatchBits<BITYPE, BITYPE_WIDTH>[numQueries]();
			}


			BatchBits<BITYPE, BITYPE_WIDTH>** next_verstatus;
			next_verstatus = new BatchBits<BITYPE, BITYPE_WIDTH>*[subgraphSize]();
			for (int a = 0; a < subgraphSize; a++) {
				const auto ret = posix_memalign(reinterpret_cast<void**>(&(next_verstatus[a])), 64, sizeof(BatchBits<BITYPE, BITYPE_WIDTH>)*numQueries);//数据对齐（分配内存首地址，对齐边界，指定分配字节大小）
				if (unlikely(ret != 0)) {
					std::cout << "unlikely" << std::endl;
					throw - 1;
				}
				new(next_verstatus[a]) BatchBits<BITYPE, BITYPE_WIDTH>[numQueries]();
			}

			std::array<Bitset*, 2> visitLists;
			for (int a = 0; a < 2; a++) {
				const auto ret = posix_memalign(reinterpret_cast<void**>(&(visitLists[a])), 64, sizeof(Bitset)*subgraphSize);//数据对齐（分配内存首地址，对齐边界，指定分配字节大小）
				if (unlikely(ret != 0)) {
					throw - 1;
				}
				new(visitLists[a]) Bitset[subgraphSize]();
			}

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
				auto per = bfsData[pos].person;
				assert(seen[per].isAllZero());
				//seen[per].setBit(pos);
				assert(!seen[per].isAllZero());
				(visitLists[0])[per].setBit(pos);
				(verstatus[per])[pos].negate();

				minPerson = std::min(minPerson, per);
#ifdef BI_DIRECTIONAl
				visitNeighbors += subgraph.retrieve(per)->size();//计算源点邻居数目和
#endif
			}

			// Initialize iteration workstate
			Bitset processQuery;//processQuery中的每一位，1表示bfs没有搜索完，0表示已经搜索完
			processQuery.negate();//全置为1

			uint32_t queriesToProcess = numQueries;//代表还剩余的bfs数
			alignas(64) uint32_t numDistDiscovered[BATCH_BITS_COUNT];//512，numDistDiscovered每一个数据代表该bfs在一轮搜索中搜索到的新的顶点数
			memset(numDistDiscovered, 0, BATCH_BITS_COUNT * sizeof(uint32_t));

			BatchDistance<bit_t, width> batchDist(numDistDiscovered);

			size_t cur_Q = 0;
			size_t curToVisitQueue = 0;
			uint32_t nextDistance = 1;//代表bfs的level

			PersonId startPerson = minPerson;

			Bitset*  toVisit = visitLists[curToVisitQueue];//visit和visitNext轮替
			Bitset*  nextToVisit = visitLists[1 - curToVisitQueue];
			//awfy::FixedSizeQueue<Quadruple<BatchBits<BITYPE, BITYPE_WIDTH>>> &Q = Q_queue0;
			bool isempty = true;
			do {
				toVisit = visitLists[curToVisitQueue];//visit和visitNext轮替
				nextToVisit = visitLists[1 - curToVisitQueue];
				isempty = true;
				size_t startTime = tschrono::now();
				assert(toVisit != nullptr);
				assert(nextToVisit != nullptr);

				/*for (int ii = 0; ii < subgraphSize; ii++) {
					std::bitset<64> bis(verstatus[ii][0].data[0]);
					//std::bitset<8> bis1(verstatus[ii].data[1]);
					std::cout << subgraph.mapInternalNodeId(ii) <<" "<< bis <<std::endl;
				}
				std::cout << std::endl;*/


				isempty = runBatchRound(subgraph, startPerson, subgraphSize, toVisit, nextToVisit,
					batchDist, processQuery, seen, nextDistance, numQueries, verstatus, next_verstatus, bfsData);

				TraceStats<500>& stats = TraceStats<500>::getStats();
				stats.traceRound(nextDistance);
				stats.addRoundDuration(nextDistance, (tschrono::now() - startTime));
				startPerson = 0;
				curToVisitQueue = 1 - curToVisitQueue;

				nextDistance++;


			} while (isempty == false);
			//free(seen);

			/*for (int ii = 0; ii < subgraphSize; ii++) {
				std::bitset<64> bis(verstatus[ii][0].data[0]);
				//std::bitset<8> bis1(verstatus[ii].data[1]);
				std::cout << subgraph.mapInternalNodeId(ii) << " " << bis << std::endl;
			}
			std::cout << std::endl;*/

			std::cout << "nextDistance = " << nextDistance << std::endl;
			for (int a = 0; a < 2; a++) {
				free(visitLists[a]);
			}
			for (int a = 0; a < subgraphSize; a++) {
				free(verstatus[a]);
			}
			for (int a = 0; a < subgraphSize; a++) {
				free(next_verstatus[a]);
			}
			//free(Q);
			//free(&Q_queue0);
#ifdef STATISTICS
			statistics.finishBatch();
#endif
		}

#ifdef SORTED_NEIGHBOR_PROCESSING
		//runBatchRound(subgraph, startPerson, subgraphSize, toVisit, nextToVisit, seen, batchDist, processQuery,Q,next_Q)
		static bool __attribute__((hot)) runBatchRound(const PersonSubgraph& subgraph, const PersonId startPerson,
			const PersonId limit, Bitset* visitList, Bitset* nextVisitList, BatchDistance<bit_t, width>& batchDist,
			const Bitset processQuery, Bitset* seen, const uint32_t nextDistance, int numQ,
			BatchBits<BITYPE, BITYPE_WIDTH>** verstatus, BatchBits<BITYPE, BITYPE_WIDTH>** next_verstatus, std::vector<BatchBFSdata>& bfsData) {

			bool is_empty = true;
			/*#ifdef DO_PREFETCH
						const int p2 = min(PREFETCH, (unsigned int)(limit - startPerson));//PREFETCH=38
						for (int a = 1; a < p2; a++) {
							__builtin_prefetch(visitList + a, 0);//visitList数据预取，不超过38，少了startperson???
							// pref=(visitList + a)->data[0];
						}
			#endif*/

			for (PersonId curPerson = startPerson; curPerson < limit; ++curPerson) {
				auto curVisit = visitList[curPerson];
				/*#ifdef DO_PREFETCH
								if (curPerson + PREFETCH < limit) {
									__builtin_prefetch(visitList + curPerson + PREFETCH, 0);
									// pref=(visitList + curPerson + PREFETCH)->data[0];
								}
				#endif*/
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
				for (int a = 0; a < numQ; a++) {
					__builtin_prefetch(&verstatus[curPerson][a], 0);
					// pref=(nextVisitList + *(friendsBounds.first+a))->data[0];
				}
#endif
				int numfriends = 0;
				int numf = 0;
				while (friendsBounds.first != friendsBounds.second) {
					const auto friends = subgraph.edgesbit[curPerson][numfriends];
#ifdef DO_PREFETCH
					//if (friendsBounds.first + PREFETCH < friendsBounds.second) {
						//__builtin_prefetch(subgraph.edgesbit[curPerson] + PREFETCH , 0);
						// pref=(nextVisitList + *(friendsBounds.first+PREFETCH))->data[0];
					//}
#endif

					const auto curFriend = visitList[*friendsBounds.first];


					for (int i = 0; i < width; i++) {
						bit_t cover = curVisit.data[i];
						for (int j = 0; j < TYPE_BITS; j++) {
							if (BitBaseOp<bit_t>::notZero(cover & BitBaseOp<bit_t>::getSetMask(j))) {
								BITYPE D;
								for (int s = 0; s < BITYPE_WIDTH; s++) {
									next_verstatus[*friendsBounds.first][i*TYPE_BITS + j].data[s] |= verstatus[curPerson][i*TYPE_BITS + j].data[s] & friends.data[s];
								}

							}
						}
						seen[*friendsBounds.first].data[i] |= cover;
					}
					++friendsBounds.first;
					++numfriends;
				}

				//std::cout << subgraph.mapInternalNodeId(curPerson) << " " << numf << "|";
			}

#ifdef BI_DIRECTIONAl
			uint32_t frontierSize = 0;
			uint64_t nextVisitNeighbors = 0;
#endif
			//size_t startchange = tschrono::now();
			for (int k = 0; k < limit; k++) {
				for (int i = 0; i < width; i++) {
					auto seenele = seen[k].data[i];
					//visitList[k].data[i] = BitBaseOp<bit_t>::zero();
					if (BitBaseOp<bit_t>::notZero(seenele)) {
						for (int j = 0; j < TYPE_BITS; j++) {
							if (BitBaseOp<bit_t>::notZero(seenele & BitBaseOp<bit_t>::getSetMask(j))) {
								int numofone = 0;
								BITYPE DD;
								for (int l = 0; l < BITYPE_WIDTH; l++) {
									DD = next_verstatus[k][i*TYPE_BITS + j].data[l] & ~verstatus[k][i*TYPE_BITS + j].data[l];

									verstatus[k][i*TYPE_BITS + j].data[l] |= next_verstatus[k][i*TYPE_BITS + j].data[l];
									numofone += BitBaseOp<BITYPE>::popCount(DD);
									next_verstatus[k][i*TYPE_BITS + j].data[l] = BitBaseOp<BITYPE>::zero();
								}
								if (numofone == 0) {
									seenele &= ~BitBaseOp<bit_t>::getSetMask(j);
									continue;
								}
								is_empty = false;
								bfsData[i*TYPE_BITS + j].totalDistances += numofone * nextDistance;
								bfsData[i*TYPE_BITS + j].totalReachable += numofone;
							}
						}
					}
					nextVisitList[k].data[i] = seenele;
					seen[k].data[i] = BitBaseOp<bit_t>::zero();
				}
			}
			//for(int a=0;a< subgraphSize)
			memset(visitList, 0, subgraph.size() * sizeof(Bitset));
			return is_empty;
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
		static void updateProcessQuery(Bitset& processQuery, const uint32_t pos, const uint32_t numDiscovered,
			BatchBFSdata& bfsData, const uint32_t distance, uint32_t& queriesToProcess) {
			auto field = pos / Bitset::TYPE_BITS_COUNT;
			auto field_bit = pos - (field*Bitset::TYPE_BITS_COUNT);

			// if(BitBaseOp<bit_t>::notZero(processQuery.data[field] & BitBaseOp<bit_t>::getSetMask(field_bit))) {
			bfsData.totalReachable += numDiscovered;//某一个源点总共搜索的顶点，除去源点
			bfsData.totalDistances += numDiscovered * distance;//某一个源点总共搜索的距离

			/*if((bfsData.componentSize-1)==bfsData.totalReachable) {
			   processQuery.data[field] = BitBaseOp<bit_t>::andNot(processQuery.data[field], BitBaseOp<bit_t>::getSetMask(field_bit));
			   queriesToProcess--;
			}
		 } else {
			assert(numDiscovered == 0);
		 }*/
		}
	};


}