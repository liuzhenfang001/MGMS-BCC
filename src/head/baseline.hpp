//Copyright (C) 2014 by Manuel Then, Moritz Kaufmann, Fernando Chirigati, Tuan-Anh Hoang-Vu, Kien Pham, Alfons Kemper, Huy T. Vo
//
//Code must not be used, distributed, without written consent by the authors
#pragma once

#include "base.hpp"
#include "statistics.hpp"
//#define COUNTBASIC
namespace Query4 {

struct BaseRunner {
	static const size_t TYPE = 4;
	static const size_t WIDTH = 1;
	static const size_t TYPE_BITS = 8;
	
	static constexpr size_t batchSize() {
		return 1;
	}
	static inline void runBatch(std::vector<BatchBFSdata>& bfsData, const PersonSubgraph& subgraph) {
		
		//LOG_PRINT("start to run" << bfsData.size());
		for (size_t i = 0; i < bfsData.size(); i++) {
			int num_neighber = 0;
			int num_stat = 0;
			//LOG_PRINT(bfsData[i].person);
			run(bfsData[i].person, subgraph, bfsData[i], num_neighber, num_stat);
		}
	}

private:
	static void run(const PersonId start, const PersonSubgraph& subgraph, BatchBFSdata& bfsData, int num_neighber, int num_stat) {
		const auto subgraphSize = subgraph.size();
		const auto sourcePer = bfsData.person;

		alignas(64) uint8_t verstatus[subgraphSize];
		memset(verstatus, 0, sizeof(verstatus));
		verstatus[sourcePer] = 1;

		alignas(64) uint8_t visit[subgraphSize];
		memset(visit, 0, sizeof(visit));
		alignas(64) uint8_t next_visit[subgraphSize];
		memset(next_visit, 0, sizeof(next_visit));
		visit[sourcePer] = 1;
		int distance = 1;
		int visit_count = 0;
		bool visit_falg = true;
		while (visit_falg) {
			visit_falg = false;
			auto cur_visit = visit_count == 0 ? visit : next_visit;
			auto cur_next_visit = visit_count == 0 ? next_visit : visit;
			for (int per = 0; per < subgraphSize; per++) {
				if (cur_visit[per] == 0)
					continue;
				const auto& curFriends = *subgraph.retrieve(per);
				auto friendsBounds = curFriends.bounds();
				while (friendsBounds.first != friendsBounds.second) {
#ifdef COUNTBASIC
					num_neighber++;
					num_stat++;
#endif // COUNTBASIC

					if (verstatus[*friendsBounds.first] == 0) {
						verstatus[*friendsBounds.first] = 1;
						bfsData.totalDistances[0] += distance;
						bfsData.totalReachable[0]++;
						cur_next_visit[*friendsBounds.first] = 1;
						visit_falg = true;
					}
					++friendsBounds.first;
				}
			}
			TraceStats<500>& stats = TraceStats<500>::getStats();
			stats.traceRound(distance);
			

#ifdef COUNTBASIC
			stats.addcount(num_neighber, num_stat);
			num_neighber = 0;
			num_stat = 0;
#endif // COUNTBASIC

			
			distance++;
			visit_count = 1 - visit_count;
			memset(cur_visit, 0, sizeof(cur_visit));
		}
	}

};

}