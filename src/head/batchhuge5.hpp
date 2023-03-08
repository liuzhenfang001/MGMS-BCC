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

#define SORTED_NEIGHBOR_PROCESSING
#define DO_PREFETCH
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
		typedef BatchBits<bit_t, width> Bitset;//Represents the BFS visited states of a vertex

		static constexpr uint64_t batchSize() {
			return BATCH_BITS_COUNT;
		}

		static void runBatch(std::vector<BatchBFSdata>& bfsData, const PersonSubgraph& subgraph
#ifdef STATISTICS
			, BatchStatistics& statistics
#endif
		) {
			const uint32_t numQueries = bfsData.size();//The number of parallel source vertices, that is the number of BFSs
			assert(numQueries > 0 && numQueries <= BATCH_BITS_COUNT);
			const auto subgraphSize = subgraph.size();


			BatchBits<BITYPE, BITYPE_WIDTH>** verstatus;
			verstatus = new BatchBits<BITYPE, BITYPE_WIDTH>*[subgraphSize]();
			for (int a = 0; a < subgraphSize; a++) {
				const auto ret = posix_memalign(reinterpret_cast<void**>(&(verstatus[a])), 64, sizeof(BatchBits<BITYPE, BITYPE_WIDTH>)*numQueries);////Data alignment (allocate first memory address, align boundary, specify allocation byte size)
				if (unlikely(ret != 0)) {
					std::cout << "unlikely" << std::endl;
					throw - 1;
				}
				new(verstatus[a]) BatchBits<BITYPE, BITYPE_WIDTH>[numQueries]();
			}

			std::array<Bitset*, 2> visitLists;
			for (int a = 0; a < 2; a++) {
				const auto ret = posix_memalign(reinterpret_cast<void**>(&(visitLists[a])), 64, sizeof(Bitset)*subgraphSize);////Data alignment (allocate first memory address, align boundary, specify allocation byte size)
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


			// Initialize active queries
			for (size_t pos = 0; pos < numQueries; pos++) {
				auto per = bfsData[pos].person;
				assert(seen[per].isAllZero());
				seen[per].setBit(pos);
				assert(!seen[per].isAllZero());
				(visitLists[0])[per].setBit(pos);
				(verstatus[per])[pos].negate();

				minPerson = std::min(minPerson, per);
			}

			// Initialize iteration workstate


			uint32_t queriesToProcess = numQueries;//
			alignas(64) uint32_t numDistDiscovered[BATCH_BITS_COUNT];//
			memset(numDistDiscovered, 0, BATCH_BITS_COUNT * sizeof(uint32_t));

			BatchDistance<bit_t, width> batchDist(numDistDiscovered);

			size_t cur_Q = 0;
			size_t curToVisitQueue = 0;
			uint32_t nextDistance = 1;//

			PersonId startPerson = minPerson;

			Bitset*  toVisit = visitLists[curToVisitQueue];//
			Bitset*  nextToVisit = visitLists[1 - curToVisitQueue];
			//awfy::FixedSizeQueue<Quadruple<BatchBits<BITYPE, BITYPE_WIDTH>>> &Q = Q_queue0;
			bool isempty = true;

			do {
				toVisit = visitLists[curToVisitQueue];//
				nextToVisit = visitLists[1 - curToVisitQueue];
				isempty = true;
				size_t startTime = tschrono::now();
				assert(toVisit != nullptr);
				assert(nextToVisit != nullptr);

				isempty = runBatchRound(subgraph, startPerson, subgraphSize, toVisit, nextToVisit, seen, nextDistance, verstatus, bfsData);

				TraceStats<BATCH_BITS_COUNT>& stats = TraceStats<BATCH_BITS_COUNT>::getStats();
				stats.traceRound(nextDistance);
				stats.addRoundDuration(nextDistance, (tschrono::now() - startTime));
				startPerson = 0;
				curToVisitQueue = 1 - curToVisitQueue;

				nextDistance++;

			} while (isempty == false);
			//free(seen);
			for (int a = 0; a < 2; a++) {
				free(visitLists[a]);
			}
			for (int a = 0; a < subgraphSize; a++) {
				free(verstatus[a]);
			}
			for (int a = 0; a < subgraphSize; a++) {
				free(next_verstatus[a]);
			}
#ifdef STATISTICS
			statistics.finishBatch();
#endif
		}

#ifdef SORTED_NEIGHBOR_PROCESSING
		//runBatchRound(subgraph, startPerson, subgraphSize, toVisit, nextToVisit, seen, batchDist, processQuery,Q,next_Q)
		static bool __attribute__((hot)) runBatchRound(const PersonSubgraph& subgraph, const PersonId startPerson,
			const PersonId limit, Bitset* visitList, Bitset* nextVisitList, Bitset* seen, const uint32_t nextDistance,
			BatchBits<BITYPE, BITYPE_WIDTH>** verstatus, std::vector<BatchBFSdata>& bfsData) {

			bool is_empty = true;

#ifdef DO_PREFETCH
			const int p2 = min(PREFETCH, (unsigned int)(limit - startPerson));//PREFETCH=38
			for (int a = 1; a < p2; a++) {
				__builtin_prefetch(visitList + a, 0);//
				// pref=(visitList + a)->data[0];
			}
#endif
			std::queue<PersonId> personqueue;
			for (PersonId curPerson = startPerson; curPerson < limit; ++curPerson) {
				auto curVisit = visitList[curPerson];
#ifdef DO_PREFETCH
				if (curPerson + PREFETCH < limit) {
					__builtin_prefetch(visitList + curPerson + PREFETCH, 0);
					// pref=(visitList + curPerson + PREFETCH)->data[0];
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

				const auto& curFriendsedges = *subgraph.retrieveedges(curPerson);
				auto friendsedgesBounds = curFriendsedges.bounds();
				const auto& curFriends = *subgraph.retrieve(curPerson);//
				auto friendsBounds = curFriends.bounds();
#ifdef DO_PREFETCH
				const int p = min(PREFETCH, (unsigned int)(friendsBounds.second - friendsBounds.first));
				for (int a = 1; a < p; a++) {
					__builtin_prefetch(nextVisitList + *(friendsBounds.first + a), 1);
					// pref=(nextVisitList + *(friendsBounds.first+a))->data[0];
				}
#endif
				while (friendsBounds.first != friendsBounds.second && friendsedgesBounds.first != friendsedgesBounds.second) {
					//
#ifdef DO_PREFETCH
					if (friendsBounds.first + PREFETCH < friendsBounds.second) {
						__builtin_prefetch(nextVisitList + *(friendsBounds.first + PREFETCH), 1);
						// pref=(nextVisitList + *(friendsBounds.first+PREFETCH))->data[0];
					}
#endif




					//bit_t cover;
					for (int i = 0; i < width; i++) {
						const bit_t cover = curVisit.data[i] & visitList[*friendsBounds.first].data[i];
						if (BitBaseOp<bit_t>::notZero(cover)) {
							for (int j = 0; j < TYPE_BITS; j++) {
								if (BitBaseOp<bit_t>::notZero(cover & BitBaseOp<bit_t>::getSetMask(j))) {
									bool insert = false;
									BatchBits<BITYPE, BITYPE_WIDTH> D;
									for (int s = 0; s < BITYPE_WIDTH; s++) {
										D.data[s] = verstatus[curPerson][i*TYPE_BITS + j].data[s] & friendsedgesBounds.first->data[s] & ~verstatus[*friendsBounds.first][i*TYPE_BITS + j].data[s];
										if (BitBaseOp<BITYPE>::notZero(D.data[s])) {
											seen[*friendsBounds.first].data[i] |= BitBaseOp<bit_t>::getSetMask(j);
											next_verstatus[*friendsBounds.first][i*TYPE_BITS + j].data[s] = D.data[s];
											is_empty = false;
											//insert = true;
											uint64_t numofone = BitBaseOp<BITYPE>::NumOfOne(D.data[s]);
											bfsData[i*TYPE_BITS + j].totalDistances += numofone * nextDistance;
											bfsData[i*TYPE_BITS + j].totalReachable += numofone;
										}
									}
								}
							}
						}
						//batchDist.updateDiscovered(cover1, i);
					}
					for (int i = 0; i < width; i++) {
						//const bit_t nocover = curVisit.data[i] & ~seen[*friendsBounds.first].data[i];
						const bit_t nocover = curVisit.data[i] & ~visitList[*friendsBounds.first].data[i];
						bit_t nocover1 = nocover;
						if (BitBaseOp<bit_t>::notZero(nocover)) {
							for (int j = 0; j < TYPE_BITS; j++) {

								if (BitBaseOp<bit_t>::notZero(nocover & BitBaseOp<bit_t>::getSetMask(j))) {
									bool insert = false;
									for (int s = 0; s < BITYPE_WIDTH; s++) {

										const BITYPE D = ((verstatus[curPerson][i*TYPE_BITS + j].data[s] & friendsedgesBounds.first->data[s]) & ~(verstatus[*friendsBounds.first][i*TYPE_BITS + j].data[s]));

										if (BitBaseOp<BITYPE>::notZero(D)) {
											insert = true;
											verstatus[*friendsBounds.first][i*TYPE_BITS + j].data[s] |= D;
											uint64_t numofone = BitBaseOp<BITYPE>::NumOfOne(D);
											bfsData[i*TYPE_BITS + j].totalDistances += numofone * nextDistance;
											bfsData[i*TYPE_BITS + j].totalReachable += numofone;
										}
									}
									if (!insert) {
										nocover1 &= ~BitBaseOp<bit_t>::getSetMask(j);
									}

								}
							}
							if (BitBaseOp<bit_t>::notZero(nocover1)) {
								nextVisitList[*friendsBounds.first].data[i] |= nocover1;
								//seen[*friendsBounds.first].data[i] |= nocover1;
								is_empty = false;
							}
						}
					}

					++friendsBounds.first;
					++friendsedgesBounds.first;
				}
				//for (int i = 0; i < width; i++) {
					//visitList[curPerson].data[i] = BitBaseOp<bit_t>::zero();//reset
				//}
			}

#ifdef BI_DIRECTIONAl
			uint32_t frontierSize = 0;
			uint64_t nextVisitNeighbors = 0;
#endif
			for (int k = 0; k < limit; k++) {
				for (int i = 0; i < width; i++) {
					visitList[k].data[i] = BitBaseOp<bit_t>::zero();
					auto seenele = seen[k].data[i];
					if (BitBaseOp<bit_t>::notZero(seenele)) {
						for (int j = 0; j < TYPE_BITS; j++) {
							if (BitBaseOp<bit_t>::notZero(seenele & BitBaseOp<bit_t>::getSetMask(j))) {
								for (int l = 0; l < BITYPE_WIDTH; l++) {
									verstatus[k][i*TYPE_BITS + j].data[l] |= next_verstatus[k][i*TYPE_BITS + j].data[l];
									next_verstatus[k][i*TYPE_BITS + j].data[l] = 0;
								}

							}
						}
						seen[k].data[i] = BitBaseOp<bit_t>::zero();
					}
				}
			}


			//batchDist.finalize();
			/*#ifdef BI_DIRECTIONAl
			return std::make_pair(frontierSize, nextVisitNeighbors);
			#else*/
			return is_empty;
			//#endif
		}

}
