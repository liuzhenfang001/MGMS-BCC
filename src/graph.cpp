/**
Copyright (C) 2023/03/08 by Zhenfang Liu, Jianxiong Ye.
Code must not be used, distributed, without written consent by the authors.
*/
#include "head/tokenizer.hpp"
#include "head/graph.hpp"
#include "head/io.hpp"
#include <random>
#include <algorithm>
#include <string.h>
using namespace std;

GraphData GraphData::loadFromPath(const std::vector<NodePair>& GraphEdges, const std::vector<float>& GraphWeight) {//reload function loadFromPath
   // Count number of persons (excluding header line)

   std::unordered_map<uint64_t,uint64_t> nodeRenaming;
   std::unordered_map<uint64_t,uint64_t> revNodeRenaming;
   uint64_t nextNodeId=0;

   auto mapExternalNodeId = [&nodeRenaming, &revNodeRenaming, &nextNodeId](uint64_t id) -> uint64_t {// Lambda experssion
      if(nodeRenaming.find(id)==nodeRenaming.end()) {
         nodeRenaming[id] = nextNodeId;//record virtual id of real id
         revNodeRenaming[nextNodeId] = id;//record real id of virtual id
         return nextNodeId++;
      } else {
         return nodeRenaming[id];//return virtual id of real id
      }
   };

   vector<NodePair> uniqueEdges(GraphEdges.size());
   size_t e=0;
   NodePair last(numeric_limits<uint64_t>::max(), numeric_limits<uint64_t>::max());
   for(const NodePair& a : GraphEdges) {
      if(last.idA!=a.idA || last.idB!=a.idB) {
         last = a;
         uniqueEdges[e++] = NodePair(mapExternalNodeId(a.idA), a.idB);
      }
   }
   uniqueEdges.resize(e);
   for(NodePair& a : uniqueEdges) {
      a.idB = mapExternalNodeId(a.idB);
   }
   return GraphData(nextNodeId, move(uniqueEdges), move(revNodeRenaming));
}

void GraphData::quick_sort(vector<NodePair>& edges,vector<float>& pro ,int l, int r)
{
	if (l < r)
	{
		//Swap(s[l], s[(l + r) / 2]); //Swap this number in the middle with the first number See Note 1
		int i = l, j = r;
		auto x = edges[l];
		auto p = pro[l];
		while (i < j)
		{
			while (i < j && (edges[j].idA > x.idA ||(edges[j].idA==x.idA && edges[j].idB>x.idB))) // Find the first number less than x from right to left
				j--;
			if (i < j) {
				edges[i] = edges[j];
				pro[i++] = pro[j];
			}

			while (i < j && (edges[i].idA < x.idA || (edges[i].idA == x.idA && edges[i].idB < x.idB))) // Find the first number greater than or equal to x from left to right
				i++;
			if (i < j) {
				edges[j] = edges[i];
				pro[j--] = pro[i];
			}
		}
		edges[i] = x;
		pro[i] = p;
		quick_sort(edges,pro, l, i - 1);
		quick_sort(edges,pro, i + 1, r);
	}
}

GraphData GraphData::sampledata(const std::vector<NodePair>& GraphEdges, const std::vector<float>& GraphWeight) {//reload function loadFromPath

	std::vector<NodePair> edges;
	std::unordered_map<uint64_t, uint64_t> nodeRenaming;
	std::unordered_map<uint64_t, uint64_t> revNodeRenaming;
	uint64_t nextNodeId = 0;

	auto mapExternalNodeId = [&nodeRenaming, &revNodeRenaming, &nextNodeId](uint64_t id) -> uint64_t {//Lambda experssion
		if (nodeRenaming.find(id) == nodeRenaming.end()) {
			nodeRenaming[id] = nextNodeId;//record virtual id of real id
			revNodeRenaming[nextNodeId] = id;//record real id of virtual id
			return nextNodeId++;
		}
		else {
			return nodeRenaming[id];//return virtual id of real id
		}
	};


	//Remove duplicates
	vector<NodePair> uniqueEdges(GraphEdges.size());
	size_t e = 0;
	NodePair last(numeric_limits<uint64_t>::max(), numeric_limits<uint64_t>::max());
	for (const NodePair& a : GraphEdges) {
		if (last.idA != a.idA || last.idB != a.idB) {
			last = a;
			uniqueEdges[e++] = NodePair(mapExternalNodeId(a.idA), a.idB);
		}
	}
	uniqueEdges.resize(e);
	for (NodePair& a : uniqueEdges) {
		a.idB = mapExternalNodeId(a.idB);
	}

	return GraphData(nextNodeId, move(uniqueEdges), move(revNodeRenaming));
}
