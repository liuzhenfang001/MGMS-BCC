/**
Copyright (C) 2023/03/08 by Zhenfang Liu, Jianxiong Ye.
Code must not be used, distributed, without written consent by the authors.
*/
// #include "./log.hpp"
#include "head/bench.hpp"
#include "head/TraceStats.hpp"
#include <iostream>
#include <stdio.h>
#define GEN_BENCH_BRANCH(X, CTYPE, WIDTH)                                                                                                                                       \
    X(batchType == sizeof(CTYPE) * 8 && batchWidth == WIDTH)                                                                                                                    \
    {                                                                                                                                                                           \
        bencher = new SpecializedBFSBenchmark<Query4::HugeBatchBfs<CTYPE, WIDTH, false>>("BatchBFS " + std::to_string(sizeof(CTYPE) * 8) + " (" + std::to_string(WIDTH) + ")"); \
        maxBatchSize = sizeof(CTYPE) * 8 * WIDTH;                                                                                                                               \
        bfsType = "MGMS-BCC";                                                                                              \
    }
//bfsType = std::to_string(sizeof(CTYPE) * 8) + "_" + std::to_string(WIDTH);



/**
Run Commend:
[executaion file] [datasets file] [number of execution times] [number of thread] [Type of algorithm] [multiples of 64] [number of source vertices] f

for example:
./runBencher ./test_queries/test.txt 1 1 MGMS-BCC 1 5 f

tips:
Sample number defaults to 1024.
*/
int main(int argc, char **argv)
{
    if (argc != 6 && argc != 7 && argc != 8)
    {
        FATAL_ERROR("Not enough parameters");
    }
    int64_t a = 0;

    std::string dataname = Queries::getdataname(std::string(argv[1]));
    // LOG_PRINT("dataset: " << dataname);
    Queries queries = Queries::loadFromFile(std::string(argv[1]));
    const int numRuns = std::stoi(std::string(argv[2]));

    // size_t bfsLimit = std::numeric_limits<uint64_t>::max();
    size_t bfsLimit = argc >= 7 ? std::stoi(std::string(argv[6])) : std::numeric_limits<uint64_t>::max();
    bool checkNumTasks = argc >= 8 && argv[7][0] == 'f' ? false : true;
    size_t numThreads = std::thread::hardware_concurrency() / 2;
    if (argc > 3)
    {
        numThreads = std::stoi(std::string(argv[3]));
    }

    // LOG_PRINT("[Main] Using " << numThreads << " threads");

    size_t maxBatchSize;
    BFSBenchmark *bencher;
    std::string bfsType;
    Query4::PARABFSRunner::setThreadNum(numThreads); // Using threads inside BFS

    if (std::string(argv[4]) == "MG-BCC")  // topk -> MG-BCC
    {
        bencher = new SpecializedBFSBenchmark<Query4::BFSRunner>("Top-KBFSRunner");
        maxBatchSize = 1;
        bfsType = "MG-BCC";
    }
    else if (std::string(argv[4]) == "MS-BCC")  // mbfs -> MS-BCC
    {
        bencher = new SpecializedBFSBenchmark<Query4::MBFSRunner<uint64_t, 1, false>>("MBFSRunner");
        maxBatchSize = 8;
        bfsType = "MS-BCC";
    }
    else if (std::string(argv[4]) == "baseline")
    {
        bencher = new SpecializedBFSBenchmark<Query4::BaseRunner>("BaseRunner");
        maxBatchSize = 1;
        bfsType = "BaseLine";
    }
    else if (std::string(argv[4]) == "GCC2")
    {
        bencher = new SpecializedBFSBenchmark<Query4::HugeBatchBfsforGCC2<uint64_t, 1, false>>("HugeBatchBfsforGCC2");
        maxBatchSize = 64;
        bfsType = "GCC2";
    }
    else if (std::string(argv[4]) == "MGMS-BCC-unOPT") // basebatch -> MGMS-BCC-unOPT
    {
        bencher = new SpecializedBFSBenchmark<Query4::HugeBatchBaseBfs<uint64_t, 1, false>>("HugeBatchBaseBfs");
        maxBatchSize = 64;
        bfsType = "MGMS-BCC-unOPT";
    }
    else if (std::string(argv[4]) == "MGMS-BCC")
    {
        //const int batchType = std::stoi(std::string(argv[4]));
        //std::cout<<"batchType="<<batchType<<endl;
        bfsType = "MGMS-BCC";
        const int batchType = 64; // batchWidth defaults to 64 (the size of a word in computer is 64).
        const int batchWidth = std::stoi(std::string(argv[5])); // the multipile of a word (64). (batchWidth defaults to 1.)
        GEN_BENCH_BRANCH(if, __m128i, 8)                        // 1024.
        GEN_BENCH_BRANCH(else if, __m128i, 4)                   // 512
        GEN_BENCH_BRANCH(else if, __m128i, 1)                   // 128
#ifdef AVX2
        GEN_BENCH_BRANCH(else if, __m256i, 2) // 512
        GEN_BENCH_BRANCH(else if, __m256i, 1) // 256
#endif
        GEN_BENCH_BRANCH(else if, uint64_t, 8)  // 512 bit
        GEN_BENCH_BRANCH(else if, uint64_t, 16) // 1024
        GEN_BENCH_BRANCH(else if, uint64_t, 1)  // 64
        GEN_BENCH_BRANCH(else if, uint32_t, 16) // 512
        GEN_BENCH_BRANCH(else if, uint32_t, 1)  // 32
        GEN_BENCH_BRANCH(else if, uint16_t, 32) // 512
        GEN_BENCH_BRANCH(else if, uint16_t, 1)  // 16
        GEN_BENCH_BRANCH(else if, uint8_t, 64)  // 512
        GEN_BENCH_BRANCH(else if, uint8_t, 1)   // 8
        else
        {
            exit(-1);
        }
    }

    // Allocate additional worker threads
    Workers workers(numThreads - 1);
    size_t numSample = SAMPLE_NUMS; //SAMPLE_NUMS defaults to 1024

    for (unsigned i = 0; i < queries.queries.size(); i++)
    {
        Query query = queries.queries[i];
        // std::cout << "[Main] Executing query " << query.dataset << endl;
        GraphFromFile grapfile = GraphFromFile();
        grapfile.init(query.dataset);
        grapfile.getSource(bfsLimit);
        if (bfsType != "MS-BFS" && bfsType != "BaseLine")
        {
            size_t startSample = tschrono::now();
            const auto personGraph = Graph<Query4::PersonId>::loadUCGFromPath(grapfile.edges, grapfile.weight);
            // LOG_PRINT("Compenents: "<<personGraph.componentSizes.size());
            // for (int i = 0;i < 10; i++) {
            // LOG_PRINT("Compenent " << i << " " << personGraph.componentSizes[i] << " nodes");
            //}
            size_t SampleTime = tschrono::now() - startSample;
            if (bfsLimit > personGraph.size())
            {
                bfsLimit = personGraph.size();
            }
            if (checkNumTasks)
            {
                auto ranges = generateTasks(bfsLimit, personGraph.size(), maxBatchSize);
                auto desiredTasks = numThreads * 3;
                if (ranges.size() < desiredTasks)
                {
                    FATAL_ERROR("[Main] Not enough tasks! #Threads=" << numThreads << ", #Tasks=" << ranges.size() << ", #DesiredTasks=" << desiredTasks << ", #maxBatchSize=" << maxBatchSize);
                }
            }
            // Run benchmark
            // std::cout << "# Benchmarking " << bencher->name << " ... " << std::endl << "# ";

            bencher->initTrace(personGraph.numVertices, personGraph.numEdges, numThreads, bfsLimit, numSample, bfsType, dataname);
            bencher->run(7, personGraph, query.reference, workers, bfsLimit, SampleTime, true);
            // std::cout << bencher->lastRuntime() << "ms ";

            // std::cout << bencher->getMinTrace() << std::endl;
        }
        else
        {

            for (int i = 0; i < SAMPLE_NUMS; i++)
            {
                size_t startSample = tschrono::now();
                const auto personGraph = Graph<Query4::PersonId>::sampleGraph(grapfile.edges, grapfile.weight);
                // LOG_PRINT("Compenents: " << personGraph.componentSizes.size());
                size_t SampleTime = tschrono::now() - startSample;
                if (bfsLimit > personGraph.size())
                {
                    bfsLimit = personGraph.size();
                }
                if (checkNumTasks)
                {
                    auto ranges = generateTasks(bfsLimit, personGraph.size(), maxBatchSize);
                    auto desiredTasks = numThreads * 3;
                    if (ranges.size() < desiredTasks)
                    {
                        FATAL_ERROR("[Main] Not enough tasks! #Threads=" << numThreads << ", #Tasks=" << ranges.size() << ", #DesiredTasks=" << desiredTasks << ", #maxBatchSize=" << maxBatchSize);
                    }
                }
                bencher->initTrace(personGraph.numVertices, personGraph.numEdges, numThreads, bfsLimit, numSample, bfsType, dataname);
                bencher->run(7, personGraph, query.reference, workers, bfsLimit, SampleTime, i == (SAMPLE_NUMS - 1));
                std::cout.flush();
            }
        }
    }
    workers.close();
    return 0;
}
