// Copyright (C) 2014 by Manuel Then, Moritz Kaufmann, Fernando Chirigati, Tuan-Anh Hoang-Vu, Kien Pham, Alfons Kemper, Huy T. Vo
//
// Code must not be used, distributed, without written consent by the authors
#pragma once

#include "log.hpp"
#include <atomic>
#include <array>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <mutex>
using namespace std;
struct Closeness
{
    size_t personId;
    // double numTraversednode;
    // double totalDistance;
    vector<double> closs;
    double closeness_value = 0;
    Closeness(size_t Id, vector<double> b)
    {
        personId = Id;
        // numTraversednode = a;
        closs = b;
        // closs= ((static_cast<double>(a) / (sizeof(BITYPE) * 8 * BITYPE_WIDTH) - 1)*(static_cast<double>(a) / (sizeof(BITYPE) * 8 * BITYPE_WIDTH) - 1)) / (static_cast<double>(totalgraphsize)*static_cast<double>(totalDistances) / (sizeof(BITYPE) * 8 * BITYPE_WIDTH));//static_cast<double>((sizeof(BITYPE) * 8 * BITYPE_WIDTH))
    }

    Closeness(){};

    void getCloseness()
    {
        for (int i = 0; i < closs.size(); i++)
        {
            closeness_value += closs[i];
        }
        closeness_value /= SAMPLE_NUMS;
    }
    /*void getCloseness(double size) {
        if (this->numTraversednode == 0) {
            this->closs = 0;
        }
        this->closs= ((static_cast<double>(numTraversednode) / (sizeof(BITYPE) * 8 * BITYPE_WIDTH))*(static_cast<double>(numTraversednode) / (sizeof(BITYPE) * 8 * BITYPE_WIDTH))) / (static_cast<double>(size)*static_cast<double>(totalDistance) / (sizeof(BITYPE) * 8 * BITYPE_WIDTH));//static_cast<double>((sizeof(BITYPE) * 8 * BITYPE_WIDTH))

    }*/
};

template <size_t maxRounds = 500>
class TraceStats
{
    size_t numVertices;
    size_t numEdges;
    size_t numTraversedEdges;
    size_t batchSize;
    size_t runnerType;
    size_t typeSize;
    size_t width;
    size_t numThreads;
    size_t maxBfs;
    std::string bfsType;
    size_t sampleTime;
    size_t BFSTime;
    std::string dataset;
    long num_neigber = 0;
    double num_stat = 0;
    size_t numR;
    std::array<std::atomic<size_t>, maxRounds> numRounds; // atomic<size_t>�����ı�������ԭ�ӻ��ģ���ԭ�Ӳ���
    std::array<std::atomic<size_t>, maxRounds> roundDurations;
    // std::array<std::array<std::atomic<size_t>, numBits>, maxRounds> roundVisitBits;
    // std::array<std::array<std::atomic<size_t>, numBits>, maxRounds> roundFriendBits;
    std::vector<Closeness> clossness;

public:
    size_t changetime;
    // TraceStats() : numRounds(), roundDurations(), roundVisitBits(), roundFriendBits(),sampleTime(0), BFSTime(0), changetime(0),clossness(0){}
    TraceStats() : numRounds(), roundDurations(), sampleTime(0), BFSTime(0), changetime(0), clossness(0)
    {
    }

public:
    static TraceStats &getStats(bool reset = false)
    {
        static TraceStats *globalStats;
        static std::mutex m;
        std::lock_guard<std::mutex> lock(m);
        if (reset && globalStats != nullptr)
        {
            globalStats = nullptr;
        }
        if (globalStats == nullptr)
        {
            globalStats = new TraceStats();
        }
        return *globalStats;
    }

    void init(size_t numVertices, size_t numEdges, size_t batchSize, size_t runnerType, size_t typeSize, size_t width, size_t numThreads, size_t maxBfs, std::string bfsType, std::string dataset)
    {
        this->numVertices = numVertices;
        this->numEdges = numEdges;
        this->numTraversedEdges = 0;
        this->batchSize = batchSize;
        this->runnerType = runnerType;
        this->typeSize = typeSize;
        this->width = width;
        this->numThreads = numThreads;
        this->maxBfs = maxBfs;
        this->bfsType = bfsType;
        this->dataset = dataset;
        this->numR = 0;
        this->clossness.reserve(maxBfs);
        // this->num_neigber = 0;
        // this->num_stat = 0;
    }

    void traceRound(size_t round)
    {
        numRounds[round]++;
        numR++;
    }

    void addcount(int neigbers, double stat)
    {
        this->num_neigber += neigbers;
        this->num_stat += stat;
    }

    int closs_size()
    {
        return clossness.size();
    }

    size_t closs_index(int i)
    {
        return clossness[i].personId;
    }

    /*void addcloss(size_t person, double numT,double totalD) {
        int i = 0;
        for (; i < clossness.size(); i++) {
            if (clossness[i].personId != person)
                continue;
            else
                break;
        }
        if (i < clossness.size()) {
            clossness[i].numTraversednode += numT;
            clossness[i].totalDistance += totalD;
        }
        else {
            //std::pair<size_t, double> a = std::make_pair(person, closs);
            clossness.push_back(Closeness(person, numT, totalD));
        }
       // if (person == clossness[0].personId)
          //  std::cout << clossness[0].personId << " " << clossness[0].totalDistance << " " << clossness[0].numTraversednode;
    }*/

    void addcloseness(size_t pId, vector<double> vec)
    {
        int i = 0;
        for (; i < clossness.size(); i++)
        {
            if (clossness[i].personId != pId)
                continue;
            else
                break;
        }
        if (i < clossness.size())
        {
            for (int j = 0; j < vec.size(); j++)
            {
                clossness[i].closs[j] += vec[j];
            }
        }
        else
        {
            clossness.push_back(Closeness(pId, vec));
        }
    }

    void printcloss(int num)
    {
        for (int i = 0; i < this->clossness.size(); i++)
        {
            this->clossness[i].getCloseness();
        }
        std::sort(clossness.begin(), clossness.end(), [](const Closeness &a, const Closeness &b) { // Lambda
            return a.closeness_value > b.closeness_value || (a.closeness_value == b.closeness_value && a.personId < b.personId);
        });
        for (int i = 0; i < num; i++)
        {
            std::cout << clossness[i].personId << "---" << clossness[i].closeness_value << " | ";
            if (i % 7 == 0)
                std::cout << std::endl;
        }
        std::cout << std::endl
                  << "###" << std::endl;
    }

    size_t getnumR()
    {
        return this->numR;
    }

    void addchange(size_t round)
    {
        changetime += round;
    }

    void addBFSTime(size_t bfstime)
    {
        BFSTime += bfstime;
    }

    void addSampleTime(size_t sampletime)
    {
        sampleTime += sampletime;
    }

    void addRoundDuration(size_t round, size_t duration)
    {
        roundDurations[round] += duration;
    }

    // void addRoundVisitBits(size_t round, size_t visitBits, size_t count) {
    // roundVisitBits[round][visitBits]+=count;
    //}

    // void addRoundFriendBits(size_t round, size_t friendBits, size_t count) {
    // roundFriendBits[round][friendBits]+=count;
    //}

    void setNumTraversedEdges(size_t c)
    {
        numTraversedEdges = c;
    }

    std::string prelude(std::string metric, size_t totalDuration)
    {
        return "[" + metric + "]\t" + std::to_string(numVertices) + "\t" + std::to_string(numEdges) + "\t" + std::to_string(numThreads) + "\t" + std::to_string(typeSize) + "\t" + std::to_string(width) + "\t" + std::to_string(batchSize) + "\t" + std::to_string(BFSTime + sampleTime) + "\t" + std::to_string(BFSTime) + "\t" + std::to_string(sampleTime) + "\t" + std::to_string(maxBfs) + "\t" + bfsType;
    }

    std::string prelude1(std::string metric, size_t totalDuration)
    {
        return "[" + metric + "]\t\t" + std::to_string(numVertices) + "\t" + std::to_string(numEdges) + "\t" + std::to_string(numTraversedEdges) + "\t" + std::to_string(batchSize) + "\t" + std::to_string(numThreads) + "\t" + "MSMBFS" + "\t" + std::to_string(typeSize) + "\t" + std::to_string(width) + "\t" + std::to_string(totalDuration) + "\t" + std::to_string(maxBfs) + "\t" + bfsType;
    }
    std::string print(size_t totalDuration)
    {
        size_t lastRound = 0;
        for (int i = 1; i < maxRounds; ++i)
        {
            if (numRounds[i] > 0)
            {
                lastRound = i;
            }
        }

        std::stringstream ss;
        // ss << "Dataset\t\tnumV\tnumE\tnumTE\tSize\tThreas\tType\tsize\twidth\tAllTime\tnumQ\tBfsType\tRounds\tnRound\tRoundTime" << std::endl;
        // ss << "Dataset\t\tnumV\tnumE\tThread\tType\tSize\tWidth\tBatch\tAllTime\tNumQ\tBfsType\tRounds\tnRound\tRoundTime" << std::endl;
        ss << prelude(dataset, totalDuration) << "\t" << lastRound << "\t" << num_neigber << "\t" << num_stat << std::endl;
        // for (int i = 1; i <= lastRound; ++i) {
        //  ss<<prelude1("TDUR",totalDuration)<<"\t"<<i<<"\t"<<numRounds[i]<<"\t"<<roundDurations[i]<<" ms"<<std::endl;
        //}

        /*#ifdef TRACE
        for (int i = 1; i < lastRound; ++i) {
           for (int j = 0; j < numBits; ++j) {
              ss<<prelude("TVBITS",totalDuration)<<"\t"<<i<<"\t"<<numRounds[i]<<"\t"<<j<<"\t"<<roundVisitBits[i][j]<<std::endl;
              ss<<prelude("TFBITS",totalDuration)<<"\t"<<i<<"\t"<<numRounds[i]<<"\t"<<j<<"\t"<<roundFriendBits[i][j]<<std::endl;
           }
        }
        #endif*/

        std::string traceStr = ss.str();
        return traceStr;
    }

    void reset()
    {
        throw;
    }
};
