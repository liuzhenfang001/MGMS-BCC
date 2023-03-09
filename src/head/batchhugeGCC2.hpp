/**
Copyright (C) 2023/03/08 by Zhenfang Liu, Jianxiong Ye.
Code must not be used, distributed, without written consent by the authors.
*/
#pragma once

#include "TraceStats.hpp"
#include "statistics.hpp"
#include "base.hpp"
#include "batchdistance.hpp"
#include <queue>
#include <array>
#include <cstring>
//#define USEPRUNR
#define SORTED_NEIGHBOR_PROCESSING
//#define DO_PREFETCH
//#define BI_DIRECTIONAl

namespace Query4 {
	// Batch part
	template<typename bit_t, uint64_t width>
	struct BatchBitsGCC {
		static const size_t TYPE_BITS_COUNT = sizeof(bit_t) * 8;

		bit_t data[width];

		BatchBitsGCC() : data() {
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

	// General HugeBatchBFS loop
	template<typename bit_t = uint64_t, uint64_t width = 1, bool detectSingle = false, uint32_t alpha = 14, uint32_t beta = 24>
	struct HugeBatchBfsforGCC2 {
		static const size_t TYPE = 3;
		static const size_t WIDTH = width;//2
		static const size_t TYPE_BITS = sizeof(bit_t) * 8;//bit_t=256
#ifdef DO_PREFETCH
		static const unsigned int PREFETCH = 2;
#endif
		static const size_t BATCH_BITS_COUNT = sizeof(bit_t)*width * 8;
		typedef BatchBitsGCC<bit_t, width> Bitset;//

		static constexpr uint64_t batchSize() {
			return BATCH_BITS_COUNT;
		}

		static void runBatch(std::vector<BatchBFSdata>& bfsData, const PersonSubgraph& subgraph
#ifdef STATISTICS
			, BatchStatistics& statistics
#endif
		) {
			const uint32_t numQueries = bfsData.size();//
			assert(numQueries > 0 && numQueries <= BATCH_BITS_COUNT);
			const auto subgraphSize = subgraph.size();


			BatchBitsGCC<BITYPE, BITYPE_WIDTH>** verstatus;
			verstatus = new BatchBitsGCC<BITYPE, BITYPE_WIDTH>*[subgraphSize]();
			for (int a = 0; a < subgraphSize; a++) {
				const auto ret = posix_memalign(reinterpret_cast<void**>(&(verstatus[a])), 64, sizeof(BatchBitsGCC<BITYPE, BITYPE_WIDTH>)*numQueries);//
				if (unlikely(ret != 0)) {
					std::cout << "unlikely" << std::endl;
					throw - 1;
				}
				new(verstatus[a]) BatchBitsGCC<BITYPE, BITYPE_WIDTH>[numQueries]();
			}


			BatchBitsGCC<BITYPE, BITYPE_WIDTH>** next_verstatus;
			next_verstatus = new BatchBitsGCC<BITYPE, BITYPE_WIDTH>*[subgraphSize]();
			for (int a = 0; a < subgraphSize; a++) {
				const auto ret = posix_memalign(reinterpret_cast<void**>(&(next_verstatus[a])), 64, sizeof(BatchBitsGCC<BITYPE, BITYPE_WIDTH>)*numQueries);//
				if (unlikely(ret != 0)) {
					std::cout << "unlikely" << std::endl;
					throw - 1;
				}
				new(next_verstatus[a]) BatchBitsGCC<BITYPE, BITYPE_WIDTH>[numQueries]();
			}

			std::array<Bitset*, 2> visitLists;
			for (int a = 0; a < 2; a++) {
				const auto ret = posix_memalign(reinterpret_cast<void**>(&(visitLists[a])), 64, sizeof(Bitset)*subgraphSize);//
				if (unlikely(ret != 0)) {
					throw - 1;
				}
				new(visitLists[a]) Bitset[subgraphSize]();
			}
#ifdef USEPRUNR
			Bitset* seen;
			{
				const auto ret = posix_memalign(reinterpret_cast<void**>(&(seen)), 64, sizeof(Bitset)*subgraphSize);//数据对齐（分配内存首地址，对齐边界，指定分配字节大小）
				if (unlikely(ret != 0)) {
					throw - 1;
				}
				new(seen) Bitset[subgraphSize]();
			}
#endif
			PersonId minPerson = std::numeric_limits<PersonId>::max();


#ifdef BI_DIRECTIONAl
			bool topDown = true;
			uint64_t unexploredEdges = subgraph.numEdges;
			uint64_t visitNeighbors = 0;
			uint32_t frontierSize = numQueries;
#endif

			// Initialize active queries
			for (size_t pos = 0; pos < numQueries; pos++) {
				auto per = bfsData[pos].person;
				(visitLists[0])[per].setBit(pos);
				(verstatus[per])[pos].negate();
				//seen[per].setBit(pos);
				minPerson = std::min(minPerson, per);
#ifdef BI_DIRECTIONAl
				visitNeighbors += subgraph.retrieve(per)->size();//
#endif
			}

			//alignas(64) uint32_t numDistDiscovered[BATCH_BITS_COUNT][SAMPLE_NUMS];//
			//memset(numDistDiscovered, 0, BATCH_BITS_COUNT * SAMPLE_NUMS * sizeof(uint32_t));

			//BatchDistance<BITYPE, BITYPE_WIDTH> batchDist[BATCH_BITS_COUNT];
			//for (int i = 0; i < BATCH_BITS_COUNT; i++) {
				//batchDist[i].numDistDiscovered = numDistDiscovered[i];
			//}

			size_t curToVisitQueue = 0;
			uint32_t nextDistance = 1;//

			PersonId startPerson = minPerson;

			Bitset*  toVisit = visitLists[curToVisitQueue];//
			Bitset*  nextToVisit = visitLists[1 - curToVisitQueue];
			bool isempty = true;
			do {
				toVisit = visitLists[curToVisitQueue];//
				nextToVisit = visitLists[1 - curToVisitQueue];
				isempty = true;
				size_t startTime = tschrono::now();
				assert(toVisit != nullptr);
				assert(nextToVisit != nullptr);

				isempty = runBatchRound(subgraph, startPerson, subgraphSize, toVisit, nextToVisit,
					nextDistance, numQueries, verstatus, next_verstatus, bfsData
#ifdef USEPRUNR
					, seen
#endif
				);
				TraceStats<500>& stats = TraceStats<500>::getStats();
				stats.traceRound(nextDistance);
				stats.addRoundDuration(nextDistance, (tschrono::now() - startTime));
				startPerson = 0;
				curToVisitQueue = 1 - curToVisitQueue;

				nextDistance++;

				//for (int i = 0; i < BATCH_BITS_COUNT; i++) {
					//batchDist[i].finalize();
				//}

			} while (isempty == false);
			//std::cout << "totalDistance = " << nextDistance << std::endl;
			for (int a = 0; a < 2; a++) {
				free(visitLists[a]);
			}
			for (int a = 0; a < subgraphSize; a++) {
				free(verstatus[a]);
			}
			for (int a = 0; a < subgraphSize; a++) {
				free(next_verstatus[a]);
			}
#ifdef USEPRUNR
			free(seen);
#endif

#ifdef STATISTICS
			statistics.finishBatch();
#endif
		}

#ifdef SORTED_NEIGHBOR_PROCESSING
		//runBatchRound(subgraph, startPerson, subgraphSize, toVisit, nextToVisit, seen, batchDist, processQuery,Q,next_Q)
		static bool __attribute__((hot)) runBatchRound(const PersonSubgraph& subgraph, const PersonId startPerson,
			const PersonId limit, Bitset* visitList, Bitset* nextVisitList,  const uint32_t nextDistance, int numQ,
			BatchBitsGCC<BITYPE, BITYPE_WIDTH>** verstatus, BatchBitsGCC<BITYPE, BITYPE_WIDTH>** next_verstatus, std::vector<BatchBFSdata>& bfsData
#ifdef USEPRUNR
			, Bitset* seen
#endif
		) {

			bool is_empty = true;
#ifdef DO_PREFETCH
			const int p2 = min(PREFETCH, (unsigned int)(limit - startPerson));//PREFETCH=38
			//for (int a = 0; a < p2; a++) {
				//__builtin_prefetch(visitList + a + startPerson, 0);//
				// pref=(visitList + a)->data[0];
			//}

			/*priorty_queue < std::pair<uint32_t, uint32_t>, std::vector<std::pair<uint32_t, uint32_t>>,
				[](std::pair<uint32_t, uint32_t> a, std::pair<uint32_t, uint32_t>, b) {
					return a.first < b.first || a.first == b.first && a.second < b.second;
				} > Blockque;*/
			for (int i = startPerson; i < startPerson + p2; i++) {
				for (int w = 0; w < WIDTH; w++) {
					for (int bit = 0; bit < TYPE_BITS; bit++) {
						if (BitBaseOp<bit_t>::notZero(visitList[i].data[w] & BitBaseOp<bit_t>::getSetMask(bit))) {
							__builtin_prefetch(verstatus[i] + w * TYPE_BITS + bit, 0);
						}
					}
				}
			}

			//for (int a = 0; a < p2; a++) {
				//__builtin_prefetch(verstatus[a], 0);
			//}

#endif

			for (PersonId curPerson = startPerson; curPerson < limit; ++curPerson) {
				auto curVisit = visitList[curPerson];

#ifdef DO_PREFETCH
				if (curPerson + PREFETCH < limit) {
					//__builtin_prefetch(visitList + curPerson + PREFETCH, 0);
					// pref=(visitList + curPerson + PREFETCH)->data[0];
					for (int w = 0; w < WIDTH; w++) {
						for (int bit = 0; bit < TYPE_BITS; bit++) {
							if (BitBaseOp<bit_t>::notZero(visitList[curPerson + PREFETCH].data[w] & BitBaseOp<bit_t>::getSetMask(bit))) {
								__builtin_prefetch(verstatus[curPerson + PREFETCH] + w * TYPE_BITS + bit, 0);
							}
						}
					}
				}
				//if (curPerson + PREFETCH < limit) {
					//__builtin_prefetch(verstatus[curPerson], 0);
				//}

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
#ifdef USEPRUNR
				auto curSeen = seen[curPerson];
#endif

				const auto& curFriends = *subgraph.retrieve(curPerson);//
				auto friendsBounds = curFriends.bounds();
#ifdef DO_PREFETCH
				//const int p = min(PREFETCH, (unsigned int)(friendsBounds.second - friendsBounds.first));

				//for (int a = 0; a < p; a++) {
					//__builtin_prefetch(nextVisitList + *(friendsBounds.first + a), 1);
					// pref=(nextVisitList + *(friendsBounds.first+a))->data[0];
				//}
				//for (int a = 0; a < p; a++) {
					//__builtin_prefetch(subgraph.edgesbit[curPerson] + a, 0);
					// pref=(nextVisitList + *(friendsBounds.first+a))->data[0];
				//}
#endif
				int numfriends = 0;
				while (friendsBounds.first != friendsBounds.second) {
					const auto friends = subgraph.edgesbit[curPerson][numfriends];
#ifdef DO_PREFETCH
					//if (friendsBounds.first + PREFETCH < friendsBounds.second) {
						//__builtin_prefetch(nextVisitList + *(friendsBounds.first+PREFETCH) , 0);
						//__builtin_prefetch(subgraph.edgesbit[curPerson] + PREFETCH, 0);
						// pref=(nextVisitList + *(friendsBounds.first+PREFETCH))->data[0];
					//}
#endif

					const auto curFriend = visitList[*friendsBounds.first];


					for (int i = 0; i < width; i++) {
#ifdef USEPRUNR
						bit_t cover = curVisit.data[i] & ~curSeen.data[i];
#else
						bit_t cover = curVisit.data[i];
#endif
						//bit_t cover = curVisit.data[i] ;
						for (int j = 0; j < TYPE_BITS; j++) {
							if (BitBaseOp<bit_t>::notZero(cover & BitBaseOp<bit_t>::getSetMask(j))) {
								//bool visit = false;
								for (int s = 0; s < BITYPE_WIDTH; s++) {
									next_verstatus[*friendsBounds.first][i*TYPE_BITS + j].data[s] |= verstatus[curPerson][i*TYPE_BITS + j].data[s] & friends.data[s];
									//if (BitBaseOp<BITYPE>::notZero(D)) {
										//next_verstatus[*friendsBounds.first][i*TYPE_BITS + j].data[s] |= D;
										//visit = true;
									//}
								}
								//if (!visit) {
									//cover = cover & ~BitBaseOp<bit_t>::getSetMask(j);
								//}

							}
						}
						if (BitBaseOp<bit_t>::notZero(cover))
							nextVisitList[*friendsBounds.first].data[i] |= cover;
					}
					++friendsBounds.first;
					++numfriends;
				}

			}

#ifdef DO_PREFETCH
			const int p3 = min(PREFETCH, (unsigned int)(limit, 0));//PREFETCH=38
			//for (int a = 0; a < p2; a++) {
				//__builtin_prefetch(visitList + a + startPerson, 0);//
				// pref=(visitList + a)->data[0];
			//}

			/*priorty_queue < std::pair<uint32_t, uint32_t>, std::vector<std::pair<uint32_t, uint32_t>>,
				[](std::pair<uint32_t, uint32_t> a, std::pair<uint32_t, uint32_t>, b) {
					return a.first < b.first || a.first == b.first && a.second < b.second;
				} > Blockque;*/
			for (int i = 0; i < p3; i++) {
				for (int w = 0; w < WIDTH; w++) {
					for (int bit = 0; bit < TYPE_BITS; bit++) {
						if (BitBaseOp<bit_t>::notZero(visitList[i].data[w] & BitBaseOp<bit_t>::getSetMask(bit))) {
							__builtin_prefetch(verstatus[i] + w * TYPE_BITS + bit, 0);
						}
					}
				}
			}

			//for (int a = 0; a < p2; a++) {
				//__builtin_prefetch(verstatus[a], 0);
			//}

#endif
			BITYPE DD;
			bool zero_one = false;
			//size_t startchange = tschrono::now();
			int numofone = 0;

			for (int k = 0; k < limit; k++) {

#ifdef DO_PREFETCH
				if (k + PREFETCH < limit) {
					//__builtin_prefetch(visitList + curPerson + PREFETCH, 0);
					// pref=(visitList + curPerson + PREFETCH)->data[0];
					for (int w = 0; w < WIDTH; w++) {
						for (int bit = 0; bit < TYPE_BITS; bit++) {
							if (BitBaseOp<bit_t>::notZero(visitList[k + PREFETCH].data[w] & BitBaseOp<bit_t>::getSetMask(bit))) {
								__builtin_prefetch(verstatus[k + PREFETCH] + w * TYPE_BITS + bit, 0);
							}
						}
					}
				}
				//if (curPerson + PREFETCH < limit) {
					//__builtin_prefetch(verstatus[curPerson], 0);
				//}

#endif
				//auto curNextVisit = nextVisitList[k];
				for (int i = 0; i < width; i++) {
					auto seenele = nextVisitList[k].data[i];
					if (BitBaseOp<bit_t>::notZero(seenele)) {
						for (int j = 0; j < TYPE_BITS; j++) {
							if (BitBaseOp<bit_t>::notZero(seenele & BitBaseOp<bit_t>::getSetMask(j))) {
								zero_one = true;
								numofone = 0;
								for (int l = 0; l < BITYPE_WIDTH; l++) {
									DD = next_verstatus[k][i*TYPE_BITS + j].data[l] & ~verstatus[k][i*TYPE_BITS + j].data[l];
									if (BitBaseOp<BITYPE>::notZero(DD)) {
										numofone += BitBaseOp<BITYPE>::popCount(DD);
										zero_one = false;
										verstatus[k][i*TYPE_BITS + j].data[l] |= next_verstatus[k][i*TYPE_BITS + j].data[l];

									}
									next_verstatus[k][i*TYPE_BITS + j].data[l] = BitBaseOp<BITYPE>::zero();//
								}
								bfsData[i*TYPE_BITS + j].totalReachable[0] += numofone;
								bfsData[i*TYPE_BITS + j].totalDistances[0] += numofone * nextDistance;
#ifdef USEPRUNR
								if (zero_one) {
									seenele &= ~BitBaseOp<bit_t>::getSetMask(j);
									seen[k].setBit(j);
									continue;
								}
#else
								if (zero_one) {
									seenele &= ~BitBaseOp<bit_t>::getSetMask(j);
									continue;
								}
#endif

								is_empty = false;
							}
						}
						nextVisitList[k].data[i] &= seenele;
					}
				}
			}
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
				__builtin_prefetch(seen + a, 0);//
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

		//updateProcessQuery(processQuery, pos, numDistDiscovered[pos], bfsData[pos], nextDistance, queriesToProcess);
		/*static void updateProcessQuery(Bitset& processQuery, const uint32_t pos, const uint32_t numDiscovered,
			BatchBFSdata& bfsData, const uint32_t distance, uint32_t& queriesToProcess) {
			auto field = pos / Bitset::TYPE_BITS_COUNT;
			auto field_bit = pos - (field*Bitset::TYPE_BITS_COUNT);

			// if(BitBaseOp<bit_t>::notZero(processQuery.data[field] & BitBaseOp<bit_t>::getSetMask(field_bit))) {
			bfsData.totalReachable += numDiscovered;//
			bfsData.totalDistances += numDiscovered * distance;//

			/*if((bfsData.componentSize-1)==bfsData.totalReachable) {
			   processQuery.data[field] = BitBaseOp<bit_t>::andNot(processQuery.data[field], BitBaseOp<bit_t>::getSetMask(field_bit));
			   queriesToProcess--;
			}
		 } else {
			assert(numDiscovered == 0);
		 }
		}*/
	};


}
