/**
Copyright (C) 2023/03/08 by Zhenfang Liu, Jianxiong Ye.
Code must not be used, distributed, without written consent by the authors.
*/
#include "tokenizer.hpp"
#include "query4.hpp"
#include "graph.hpp"
//#include "query4.hpp"
#include "noqueue.hpp"
#include "io.hpp"
#include "worker.hpp"
#include "TraceStats.hpp"
#include "parabfs.hpp"

struct Query {
   uint64_t numNodes;
   std::string dataset;
   std::string reference;
};

struct Queries {
   std::vector<Query> queries;
   static std::string getdataname(const std::string& dataset) {
	   char path_str[100];
	   strcpy(path_str, dataset.c_str());
	   char *p;
	   p = strtok(path_str, "/");
	   p = strtok(NULL, ".");
	   std::string dataname = p;
	   return dataname;
   }

   static Queries loadFromFile(const std::string& queryFile) {
	  std::string re_file = "../" + queryFile;
      io::MmapedFile file(re_file, O_RDONLY);

      Queries queries;

      Tokenizer tokenizer(file.mapping, file.size);
      // Read queries
      while(!tokenizer.isFinished()) {
         Query query;
         query.numNodes = tokenizer.readId(' ');//TXT file first data "id"
         query.dataset = tokenizer.readStr(' ');//Obtain the csv file name
         query.reference = tokenizer.readStr('\n');//Get everything in the file (top-k reference of closeness centrality)
         if(query.reference.size()==0) {
            FATAL_ERROR("[Queries] Could not load reference result, is it specified?");
         }
         queries.queries.push_back(query);
      }

      //LOG_PRINT("[Queries] Loaded "<< queries.queries.size()<< " queries.");
      return queries;
   }
};

struct BFSBenchmark {
   const std::string name;
   std::vector<uint32_t> runtimes;

   BFSBenchmark(std::string name)
      : name(name)
   { }
   virtual ~BFSBenchmark() { }

   uint32_t lastRuntime() const {
      return runtimes.back();
   }

   uint32_t avgRuntime() const {
      uint32_t sum=0;
      for(auto t : runtimes) {
         sum += t;
      }
      return sum / runtimes.size();
   }

   virtual void run(const uint32_t k, const Query4::PersonSubgraph& subgraph, const string& referenceResult, Workers& workers, uint64_t maxBfs, size_t sampletime,bool finish) = 0;

   virtual size_t batchSize() = 0;

   virtual void initTrace(size_t numVertices, size_t numEdges, size_t numThreads, size_t maxBfs, std::string bfsType, std::string dataname) = 0;

   virtual std::string getMinTrace() = 0;
};

template<typename BFSRunnerT>
struct SpecializedBFSBenchmark : public BFSBenchmark {
   #ifdef STATISTICS
   Query4::BatchStatistics statistics;
   #endif

   std::vector<std::string> traces;

   typedef TraceStats<500> RunnerTraceStats;

   SpecializedBFSBenchmark(std::string name)
      : BFSBenchmark(name)
   { }
   virtual void run(const uint32_t k, const Query4::PersonSubgraph& subgraph, const string& referenceResult, Workers& workers, uint64_t maxBfs,size_t sampletime,bool finish) override {
      uint64_t runtime;
      std::string result = runBFS<BFSRunnerT>(k, subgraph, workers, maxBfs, runtime
      #ifdef STATISTICS
         ,statistics
         #endif
         );
      //if(maxBfs == std::numeric_limits<uint64_t>::max() && result != referenceResult) {
         //cout<<endl;
         //FATAL_ERROR("[Query] Wrong result, expected ["<<referenceResult<<"], got ["<<result<<"]");
      //}
      runtimes.push_back(runtime);
	  //std::cout << result << std::endl;
      RunnerTraceStats& stats = RunnerTraceStats::getStats();
	  //std::cout << sampletime << " " << runtime << std::endl;
	  stats.addSampleTime(sampletime);
	  stats.addBFSTime(runtime);
	  //std::cout << stats.getnumR() << std::endl;
	  if(finish){
		  //std::cout<<stats.closs_size()<<std::endl;
		  stats.printcloss(maxBfs);
		  FILE *f = NULL;
		  f = fopen("./result.txt", "a");
		  fprintf(f, (stats.print(runtime)).c_str());
		  fclose(f);
		  traces.push_back(stats.print(runtime));
	  }
   }

   virtual size_t batchSize() {
      return BFSRunnerT::batchSize();
   }

   virtual void initTrace(size_t numVertices, size_t numEdges, size_t numThreads, size_t maxBfs, std::string bfsType, std::string dataname) {
      RunnerTraceStats& stats = RunnerTraceStats::getStats();
      // TODO: Account for overriding of batch size by env variale
      stats.init(numVertices, numEdges, batchSize(), BFSRunnerT::TYPE, BFSRunnerT::TYPE_BITS, BFSRunnerT::WIDTH, numThreads, maxBfs, bfsType, dataname);
   }

   virtual std::string getMinTrace() {
      size_t minRuntimeIx=0;
      for (int i = 0; i < runtimes.size(); ++i) {
         if(runtimes[i]<runtimes[minRuntimeIx]) {
            minRuntimeIx=i;
         }
      }

      return traces[minRuntimeIx];
   }
};
