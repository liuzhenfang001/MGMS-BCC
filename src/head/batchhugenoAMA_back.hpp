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
// #define USEPRUNR
#ifndef USE_BATCHDIST
#define USE_BATCHDIST
#endif // !USE_BATCH

// #define USE_BATCHDIST
#define SORTED_NEIGHBOR_PROCESSING
// #define DO_PREFETCH
// #define BI_DIRECTIONAl

namespace Query4
{
    // Batch part
    template <typename bit_t, uint64_t width>
    struct BatchBitsAMA
    {
        static const size_t TYPE_BITS_COUNT = sizeof(bit_t) * 8;

        bit_t data[width];

        BatchBitsAMA() : data()
        {
        }

#ifdef DEBUG
        bool isAllZero() const
        {
            for (unsigned i = 0; i < width; ++i)
            {
                if (BitBaseOp<bit_t>::notZero(data[i]))
                {
                    return false;
                }
            }

            return true;
        }
#endif

#ifdef STATISTICS
        size_t count() const
        {
            size_t sum = 0;
            for (unsigned i = 0; i < width; ++i)
            {
                sum += BitBaseOp<bit_t>::popCount(data[i]);
            }
            return sum;
        }
#endif

        void negate()
        {
            for (unsigned i = 0; i < width; ++i)
            {
                data[i] = ~data[i];
            }
        }

        void setBit(const size_t ix)
        {
            auto field = ix / TYPE_BITS_COUNT;
            auto field_bit = ix - (field * TYPE_BITS_COUNT);
            data[field] |= BitBaseOp<bit_t>::getSetMask(field_bit);
        }
    };

    // General HugeBatchBFS loop
    template <typename bit_t = uint64_t, uint64_t width = 1, bool detectSingle = false, uint32_t alpha = 14, uint32_t beta = 24>
    struct simpleHugeBatchBfs
    {
        static const size_t TYPE = 3;
        static const size_t WIDTH = width;                 // 2
        static const size_t TYPE_BITS = sizeof(bit_t) * 8; // bit_t=256
#ifdef DO_PREFETCH
        static const unsigned int PREFETCH = 2;
#endif
        static const size_t BATCH_BITS_COUNT = sizeof(bit_t) * width * 8;
        typedef BatchBitsAMA<bit_t, width> Bitset; //

        static constexpr uint64_t batchSize()
        {
            return BATCH_BITS_COUNT;
        }

        static void runBatch(std::vector<BatchBFSdata> &bfsData, const PersonSubgraph &subgraph
#ifdef STATISTICS
                             ,
                             BatchStatistics &statistics
#endif
        )
        {
            const uint32_t numQueries = bfsData.size(); //
            assert(numQueries > 0 && numQueries <= BATCH_BITS_COUNT);
            const auto subgraphSize = subgraph.size();

            BatchBitsAMA<BITYPE, BITYPE_WIDTH> **verstatus;
            verstatus = new BatchBitsAMA<BITYPE, BITYPE_WIDTH> *[subgraphSize]();
            for (int a = 0; a < subgraphSize; a++)
            {
                const auto ret = posix_memalign(reinterpret_cast<void **>(&(verstatus[a])), 64, sizeof(BatchBitsAMA<BITYPE, BITYPE_WIDTH>) * numQueries); //
                if (unlikely(ret != 0))
                {
                    std::cout << "unlikely" << std::endl;
                    throw -1;
                }
                new (verstatus[a]) BatchBitsAMA<BITYPE, BITYPE_WIDTH>[ numQueries ]();
            }

            BatchBitsAMA<BITYPE, BITYPE_WIDTH> **next_verstatus;
            next_verstatus = new BatchBitsAMA<BITYPE, BITYPE_WIDTH> *[subgraphSize]();
            for (int a = 0; a < subgraphSize; a++)
            {
                const auto ret = posix_memalign(reinterpret_cast<void **>(&(next_verstatus[a])), 64, sizeof(BatchBitsAMA<BITYPE, BITYPE_WIDTH>) * numQueries); //
                if (unlikely(ret != 0))
                {
                    std::cout << "unlikely" << std::endl;
                    throw -1;
                }
                new (next_verstatus[a]) BatchBitsAMA<BITYPE, BITYPE_WIDTH>[ numQueries ]();
            }

            std::array<Bitset *, 2> visitLists;
            for (int a = 0; a < 2; a++)
            {
                const auto ret = posix_memalign(reinterpret_cast<void **>(&(visitLists[a])), 64, sizeof(Bitset) * subgraphSize); //
                if (unlikely(ret != 0))
                {
                    throw -1;
                }
                new (visitLists[a]) Bitset[subgraphSize]();
            }
#ifdef USEPRUNR
            Bitset *seen;
            {
                const auto ret = posix_memalign(reinterpret_cast<void **>(&(seen)), 64, sizeof(Bitset) * subgraphSize); // ���ݶ��루�����ڴ��׵�ַ������߽磬ָ�������ֽڴ�С��
                if (unlikely(ret != 0))
                {
                    throw -1;
                }
                new (seen) Bitset[subgraphSize]();
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
            for (size_t pos = 0; pos < numQueries; pos++)
            {
                auto per = bfsData[pos].person;
                (visitLists[0])[per].setBit(pos);
                (verstatus[per])[pos].negate();
                // seen[per].setBit(pos);
                minPerson = std::min(minPerson, per);
#ifdef BI_DIRECTIONAl
                visitNeighbors += subgraph.retrieve(per)->size(); //
#endif
            }

            alignas(64) uint32_t numDistDiscovered[BATCH_BITS_COUNT][SAMPLE_NUMS]; //
            memset(numDistDiscovered, 0, BATCH_BITS_COUNT * SAMPLE_NUMS * sizeof(uint32_t));

            BatchDistance<BITYPE, BITYPE_WIDTH> batchDist[BATCH_BITS_COUNT];
            for (int i = 0; i < BATCH_BITS_COUNT; i++)
            {
                batchDist[i].numDistDiscovered = numDistDiscovered[i];
            }

            size_t curToVisitQueue = 0;
            uint32_t nextDistance = 1; //

            PersonId startPerson = minPerson;

            Bitset *toVisit = visitLists[curToVisitQueue]; //
            Bitset *nextToVisit = visitLists[1 - curToVisitQueue];
            bool isempty = true;
            do
            {
                toVisit = visitLists[curToVisitQueue]; //
                nextToVisit = visitLists[1 - curToVisitQueue];
                isempty = true;
                size_t startTime = tschrono::now();
                assert(toVisit != nullptr);
                assert(nextToVisit != nullptr);

                isempty = runBatchRound(subgraph, startPerson, subgraphSize, toVisit, nextToVisit, batchDist,
                                        nextDistance, numQueries, verstatus, next_verstatus, bfsData
#ifdef USEPRUNR
                                        ,
                                        seen
#endif
                );
#ifdef USE_BATCHDIST
                for (int a = 0; a < BATCH_BITS_COUNT; a++)
                {
                    for (int l = 0; l < SAMPLE_NUMS; l++)
                    {
                        batchDist[a].finalize();
                        bfsData[a].totalReachable[l] += numDistDiscovered[a][l];
                        bfsData[a].totalDistances[l] += numDistDiscovered[a][l] * nextDistance;
                    }
                }
                memset(numDistDiscovered, 0, BATCH_BITS_COUNT * SAMPLE_NUMS * sizeof(uint32_t));
#endif // USE_BATCHDIST

                TraceStats<500> &stats = TraceStats<500>::getStats();
                stats.traceRound(nextDistance);
                stats.addRoundDuration(nextDistance, (tschrono::now() - startTime));
                startPerson = 0;
                curToVisitQueue = 1 - curToVisitQueue;

                nextDistance++;

                // for (int i = 0; i < BATCH_BITS_COUNT; i++) {
                // batchDist[i].finalize();
                //}

            } while (isempty == false);
            // std::cout << "totalDistance = " << nextDistance << std::endl;
            for (int a = 0; a < 2; a++)
            {
                free(visitLists[a]);
            }
            for (int a = 0; a < subgraphSize; a++)
            {
                free(verstatus[a]);
            }
            for (int a = 0; a < subgraphSize; a++)
            {
                free(next_verstatus[a]);
            }
#ifdef USEPRUNR
            free(seen);
#endif

#ifdef STATISTICS
            statistics.finishBatch();
#endif
        }

        // runBatchRound(subgraph, startPerson, subgraphSize, toVisit, nextToVisit, seen, batchDist, processQuery,Q,next_Q)
        static bool __attribute__((hot)) runBatchRound(const PersonSubgraph &subgraph, const PersonId startPerson,
                                                       const PersonId limit, Bitset *visitList, Bitset *nextVisitList, BatchDistance<BITYPE, BITYPE_WIDTH> *batchDist, const uint32_t nextDistance, int numQ,
                                                       BatchBitsAMA<BITYPE, BITYPE_WIDTH> **verstatus, BatchBitsAMA<BITYPE, BITYPE_WIDTH> **next_verstatus, std::vector<BatchBFSdata> &bfsData
#ifdef USEPRUNR
                                                       ,
                                                       Bitset *seen
#endif
        )
        {

            bool is_empty = true;
#ifdef DO_PREFETCH
            const int p2 = min(PREFETCH, (unsigned int)(limit - startPerson)); // PREFETCH=38
            // for (int a = 0; a < p2; a++) {
            //__builtin_prefetch(visitList + a + startPerson, 0);//
            //  pref=(visitList + a)->data[0];
            //}

            /*priorty_queue < std::pair<uint32_t, uint32_t>, std::vector<std::pair<uint32_t, uint32_t>>,
                [](std::pair<uint32_t, uint32_t> a, std::pair<uint32_t, uint32_t>, b) {
                    return a.first < b.first || a.first == b.first && a.second < b.second;
                } > Blockque;*/
            for (int i = startPerson; i < startPerson + p2; i++)
            {
                for (int w = 0; w < WIDTH; w++)
                {
                    for (int bit = 0; bit < TYPE_BITS; bit++)
                    {
                        if (BitBaseOp<bit_t>::notZero(visitList[i].data[w] & BitBaseOp<bit_t>::getSetMask(bit)))
                        {
                            __builtin_prefetch(verstatus[i] + w * TYPE_BITS + bit, 0);
                        }
                    }
                }
            }

            // for (int a = 0; a < p2; a++) {
            //__builtin_prefetch(verstatus[a], 0);
            //}

#endif

            for (PersonId curPerson = startPerson; curPerson < limit; ++curPerson)
            {
                auto curVisit = visitList[curPerson];

#ifdef DO_PREFETCH
                if (curPerson + PREFETCH < limit)
                {
                    //__builtin_prefetch(visitList + curPerson + PREFETCH, 0);
                    // pref=(visitList + curPerson + PREFETCH)->data[0];
                    for (int w = 0; w < WIDTH; w++)
                    {
                        for (int bit = 0; bit < TYPE_BITS; bit++)
                        {
                            if (BitBaseOp<bit_t>::notZero(visitList[curPerson + PREFETCH].data[w] & BitBaseOp<bit_t>::getSetMask(bit)))
                            {
                                __builtin_prefetch(verstatus[curPerson + PREFETCH] + w * TYPE_BITS + bit, 0);
                            }
                        }
                    }
                }
                // if (curPerson + PREFETCH < limit) {
                //__builtin_prefetch(verstatus[curPerson], 0);
                //}

#endif
                bool zero = true;
                for (int i = 0; i < width; i++)
                {
                    if (BitBaseOp<bit_t>::notZero(curVisit.data[i]))
                    {
                        zero = false;
                        break;
                    }
                }
                if (zero)
                {
                    continue;
                }
#ifdef USEPRUNR
                auto curSeen = seen[curPerson];
#endif

                const auto &curFriends = *subgraph.retrieve(curPerson); //
                auto friendsBounds = curFriends.bounds();
#ifdef DO_PREFETCH
                // const int p = min(PREFETCH, (unsigned int)(friendsBounds.second - friendsBounds.first));

                // for (int a = 0; a < p; a++) {
                //__builtin_prefetch(nextVisitList + *(friendsBounds.first + a), 1);
                //  pref=(nextVisitList + *(friendsBounds.first+a))->data[0];
                //}
                // for (int a = 0; a < p; a++) {
                //__builtin_prefetch(subgraph.edgesbit[curPerson] + a, 0);
                // pref=(nextVisitList + *(friendsBounds.first+a))->data[0];
                //}
#endif
                int numfriends = 0;
                BITYPE DD;
                while (friendsBounds.first != friendsBounds.second)
                {
                    const auto friends = subgraph.edgesbit[curPerson][numfriends];
#ifdef DO_PREFETCH
                    // if (friendsBounds.first + PREFETCH < friendsBounds.second) {
                    //__builtin_prefetch(nextVisitList + *(friendsBounds.first+PREFETCH) , 0);
                    //__builtin_prefetch(subgraph.edgesbit[curPerson] + PREFETCH, 0);
                    //  pref=(nextVisitList + *(friendsBounds.first+PREFETCH))->data[0];
                    //}
#endif

                    const auto curFriend = visitList[*friendsBounds.first];

                    for (int i = 0; i < width; i++)
                    {
#ifdef USEPRUNR
                        bit_t cover = curVisit.data[i] & ~curSeen.data[i];
#else
                        bit_t cover = curVisit.data[i];
#endif
                        // bit_t cover = curVisit.data[i] ;
                        for (int j = 0; j < TYPE_BITS; j++)
                        {
                            if (BitBaseOp<bit_t>::notZero(cover & BitBaseOp<bit_t>::getSetMask(j)))
                            {
                                // int numofone = 0;
                                bool visit = false;
                                for (int s = 0; s < BITYPE_WIDTH; s++)
                                {
                                    DD = verstatus[curPerson][i * TYPE_BITS + j].data[s] & friends.data[s] & ~verstatus[*friendsBounds.first][i * TYPE_BITS + j].data[s];

                                    if (BitBaseOp<BITYPE>::notZero(DD))
                                    {
                                        visit = true;
                                        // numofone += BitBaseOp<BITYPE>::popCount(DD);
                                        batchDist[i * TYPE_BITS + j].updateDiscovered(DD, s);
                                        next_verstatus[*friendsBounds.first][i * TYPE_BITS + j].data[s] |= DD;
                                    }
                                    // next_verstatus[*friendsBounds.first][i*TYPE_BITS + j].data[s] |= verstatus[curPerson][i*TYPE_BITS + j].data[s] & friends.data[s];

                                    // if (BitBaseOp<BITYPE>::notZero(D)) {
                                    // next_verstatus[*friendsBounds.first][i*TYPE_BITS + j].data[s] |= D;
                                    // visit = true;
                                    //}
                                }
                                if (!visit)
                                {
                                    cover = cover & ~BitBaseOp<bit_t>::getSetMask(j);
                                }
                            }
                        }
                        if (BitBaseOp<bit_t>::notZero(cover))
                        {
                            is_empty = false;
                            nextVisitList[*friendsBounds.first].data[i] |= cover;
                        }
                    }
                    ++friendsBounds.first;
                    ++numfriends;
                }
            }

#ifdef BI_DIRECTIONAl
            uint32_t frontierSize = 0;
            uint64_t nextVisitNeighbors = 0;
#endif

#ifdef DO_PREFETCH
            const int p3 = min(PREFETCH, (unsigned int)(limit, 0)); // PREFETCH=38
            // for (int a = 0; a < p2; a++) {
            //__builtin_prefetch(visitList + a + startPerson, 0);//
            //  pref=(visitList + a)->data[0];
            //}

            /*priorty_queue < std::pair<uint32_t, uint32_t>, std::vector<std::pair<uint32_t, uint32_t>>,
                [](std::pair<uint32_t, uint32_t> a, std::pair<uint32_t, uint32_t>, b) {
                    return a.first < b.first || a.first == b.first && a.second < b.second;
                } > Blockque;*/
            for (int i = 0; i < p3; i++)
            {
                for (int w = 0; w < WIDTH; w++)
                {
                    for (int bit = 0; bit < TYPE_BITS; bit++)
                    {
                        if (BitBaseOp<bit_t>::notZero(visitList[i].data[w] & BitBaseOp<bit_t>::getSetMask(bit)))
                        {
                            __builtin_prefetch(verstatus[i] + w * TYPE_BITS + bit, 0);
                        }
                    }
                }
            }

            // for (int a = 0; a < p2; a++) {
            //__builtin_prefetch(verstatus[a], 0);
            //}

#endif

            // size_t startchange = tschrono::now();
            int numofone = 0;

#ifndef USE_BATCHDIST
            for (int k = 0; k < limit; k++)
            {

#ifdef DO_PREFETCH
                if (k + PREFETCH < limit)
                {
                    //__builtin_prefetch(visitList + curPerson + PREFETCH, 0);
                    // pref=(visitList + curPerson + PREFETCH)->data[0];
                    for (int w = 0; w < WIDTH; w++)
                    {
                        for (int bit = 0; bit < TYPE_BITS; bit++)
                        {
                            if (BitBaseOp<bit_t>::notZero(visitList[k + PREFETCH].data[w] & BitBaseOp<bit_t>::getSetMask(bit)))
                            {
                                __builtin_prefetch(verstatus[k + PREFETCH] + w * TYPE_BITS + bit, 0);
                            }
                        }
                    }
                }
                // if (curPerson + PREFETCH < limit) {
                //__builtin_prefetch(verstatus[curPerson], 0);
                //}

#endif
                // auto curNextVisit = nextVisitList[k];
                for (int i = 0; i < width; i++)
                {
                    auto seenele = nextVisitList[k].data[i];
                    // visitList[k].data[i] = BitBaseOp<bit_t>::zero();
                    if (BitBaseOp<bit_t>::notZero(seenele))
                    {
                        for (int j = 0; j < TYPE_BITS; j++)
                        {
                            if (BitBaseOp<bit_t>::notZero(seenele & BitBaseOp<bit_t>::getSetMask(j)))
                            {
                                zero_one = true;
                                numofone = 0;
                                for (int l = 0; l < BITYPE_WIDTH; l++)
                                {
                                    DD = next_verstatus[k][i * TYPE_BITS + j].data[l] & ~verstatus[k][i * TYPE_BITS + j].data[l];

                                    if (BitBaseOp<BITYPE>::notZero(DD))
                                    {
                                        for (int p = 0; p < sizeof(BITYPE) * 8; p++)
                                        {
                                            numofone = 0;
                                            numofone = BitBaseOp<BITYPE>::notZero(DD & BitBaseOp<BITYPE>::getSetMask(p)) ? 1 : 0;
                                            bfsData[i * TYPE_BITS + j].totalReachable[l * sizeof(BITYPE) * 8 + p] += numofone;

                                            bfsData[i * TYPE_BITS + j].totalDistances[l * sizeof(BITYPE) * 8 + p] += numofone * nextDistance;
                                        }
#ifdef USE_BATCHDIST
                                        batchDist[i * TYPE_BITS + j].updateDiscovered(DD, l);
#endif // USE_BATCHDIST
                                        zero_one = false;
                                        verstatus[k][i * TYPE_BITS + j].data[l] |= next_verstatus[k][i * TYPE_BITS + j].data[l];
                                        // numofone += BitBaseOp<BITYPE>::popCount(DD);
                                    }
                                    // batchDist[i*TYPE_BITS + j].updateDiscovered(DD,l);//

                                    next_verstatus[k][i * TYPE_BITS + j].data[l] = BitBaseOp<BITYPE>::zero(); //
                                }
                                // cout << numofone << " ";
#ifdef USEPRUNR
                                if (zero_one)
                                {
                                    seenele &= ~BitBaseOp<bit_t>::getSetMask(j);
                                    seen[k].setBit(j);
                                    continue;
                                }
#else
                                if (zero_one)
                                {
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
#else
            for (int k = 0; k < limit; k++)
            {

#ifdef DO_PREFETCH
                if (k + PREFETCH < limit)
                {
                    //__builtin_prefetch(visitList + curPerson + PREFETCH, 0);
                    // pref=(visitList + curPerson + PREFETCH)->data[0];
                    for (int w = 0; w < WIDTH; w++)
                    {
                        for (int bit = 0; bit < TYPE_BITS; bit++)
                        {
                            if (BitBaseOp<bit_t>::notZero(visitList[k + PREFETCH].data[w] & BitBaseOp<bit_t>::getSetMask(bit)))
                            {
                                __builtin_prefetch(verstatus[k + PREFETCH] + w * TYPE_BITS + bit, 0);
                            }
                        }
                    }
                }
                // if (curPerson + PREFETCH < limit) {
                //__builtin_prefetch(verstatus[curPerson], 0);
                //}

#endif
                // auto curNextVisit = nextVisitList[k];
                for (int i = 0; i < width; i++)
                {
                    auto seenele = nextVisitList[k].data[i];
                    if (BitBaseOp<bit_t>::notZero(seenele))
                    {
                        for (int j = 0; j < TYPE_BITS; j++)
                        {
                            if (BitBaseOp<bit_t>::notZero(seenele & BitBaseOp<bit_t>::getSetMask(j)))
                            {
                                for (int l = 0; l < BITYPE_WIDTH; l++)
                                {
                                    verstatus[k][i * TYPE_BITS + j].data[l] |= next_verstatus[k][i * TYPE_BITS + j].data[l];
                                    // batchDist[i*TYPE_BITS + j].updateDiscovered(DD,l);//
                                    next_verstatus[k][i * TYPE_BITS + j].data[l] = BitBaseOp<BITYPE>::zero(); //
                                }
                            }
                        }
                    }
                }
            }
#endif // USE_BATCHDIST
            memset(visitList, 0, subgraph.size() * sizeof(Bitset));
            return is_empty;
        }
    };
}
