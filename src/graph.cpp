//Copyright (C) 2014 by Manuel Then, Moritz Kaufmann, Fernando Chirigati, Tuan-Anh Hoang-Vu, Kien Pham, Alfons Kemper, Huy T. Vo
//
//Code must not be used, distributed, without written consent by the authors
#include "head/tokenizer.hpp"
#include "head/graph.hpp"
#include "head/io.hpp"
#include <random>
#include <algorithm>
#include <string.h>
using namespace std;

GraphData GraphData::loadFromPath(const std::vector<NodePair>& GraphEdges, const std::vector<float>& GraphWeight) {//重载loadFromPath函数
   // Count number of persons (excluding header line)
	
   std::unordered_map<uint64_t,uint64_t> nodeRenaming;
   std::unordered_map<uint64_t,uint64_t> revNodeRenaming;
   uint64_t nextNodeId=0;

   auto mapExternalNodeId = [&nodeRenaming, &revNodeRenaming, &nextNodeId](uint64_t id) -> uint64_t {//Lambda表达式
      if(nodeRenaming.find(id)==nodeRenaming.end()) {
         nodeRenaming[id] = nextNodeId;//记录id对应的伪id
         revNodeRenaming[nextNodeId] = id;//记录伪id对应的真实id
         return nextNodeId++;
      } else {
         return nodeRenaming[id];//返回id对应的伪id
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


/*GraphData GraphData::loadFromPath(const std::string& edgesFile, bool isUCG) {//重载loadFromPath函数
   // Count number of persons (excluding header line)
	std::string data_path = "../" + edgesFile;
	size_t numEdges = io::fileLines(data_path) - 1;

	const char *sep = "/";
	char path_str[100];
	strcpy(path_str, edgesFile.c_str());
	char *p;
	p = strtok(path_str, sep);
	p = strtok(NULL, sep);
	p = strtok(NULL, sep);
	string filename = p;
	std::string sample_data = "../test_queries/data/sample_" + filename;

	int fd = open(sample_data.c_str(), O_RDONLY);
	if (fd < 0) {
		FILE *fp;
		if ((fp = fopen(sample_data.c_str(), "w")) != NULL) {
			std::random_device rd;
			std::mt19937 g(rd());
			std::uniform_int_distribution<> dis(1, 100);
			for (int i = 0; i < numEdges; i++) {
				std::cout << "succ" << std::endl;
				fprintf(fp, "%1.2f\n", float(dis(g)) / 100);
			}
		}
		fclose(fp);
		//新建一个sample_文件，每行写上概率。
	}

	vector<float> pro;
	pro.reserve(numEdges);
	io::MmapedFile sample_file(sample_data, O_RDONLY);
	Tokenizer sample_tokenizer(sample_file.mapping, sample_file.size);
	while (!sample_tokenizer.isFinished()) {
		assert(pro.size() < numEdges * 2);
		float tem_pro = atof((sample_tokenizer.readStr('\n')).c_str());
		pro.push_back(tem_pro);
	}

	if (isUCG)
		LOG_PRINT("[LOADING] Number of edges: " << numEdges);
	// Load edges data from file
	if (isUCG)
		LOG_PRINT("[LOADING] Loading edges from file: " << data_path);
	io::MmapedFile file(data_path, O_RDONLY);

	Tokenizer tokenizer(file.mapping, file.size);
	tokenizer.skipLine(); // Skip header line 

	vector<NodePair> edges;
	edges.reserve(numEdges);

	std::unordered_map<uint64_t, uint64_t> nodeRenaming;
	std::unordered_map<uint64_t, uint64_t> revNodeRenaming;
	uint64_t nextNodeId = 0;

	auto mapExternalNodeId = [&nodeRenaming, &revNodeRenaming, &nextNodeId](uint64_t id) -> uint64_t {//Lambda表达式
		if (nodeRenaming.find(id) == nodeRenaming.end()) {
			nodeRenaming[id] = nextNodeId;//记录id对应的伪id
			revNodeRenaming[nextNodeId] = id;//记录伪id对应的真实id
			return nextNodeId++;
		}
		else {
			return nodeRenaming[id];//返回id对应的伪id
		}
	};

	while (!tokenizer.isFinished()) {
		assert(edges.size() < numEdges * 2);
		NodePair pair;
		pair.idA = tokenizer.readId('|');
		pair.idB = tokenizer.readId('\n');
		if (pair.idA == pair.idB) {
			continue; //No self-edges
		}

		//Add undirected（正向与反向)
		edges.push_back(pair);
		if (!isUCG)
			edges.push_back(NodePair(pair.idB, pair.idA));
	}
	//LOG_PRINT("size() of edges = " << edges.size());

	// Reading edges
	//LOG_PRINT("[LOADING] Read edges");
	//std::sort(edges.begin(), edges.end(), [](const NodePair& a, const NodePair& b) {//Lambda
	   //return a.idA<b.idA||(a.idA==b.idA && a.idB<b.idB);
	//});

	quick_sort(edges, pro, 0, numEdges - 1);

	//for (int i = 0; i < edges.size(); i++) {
	  // std::cout << edges[i].idA << " " << edges[i].idB<<" | ";
	   //if (i % 20 == 0)
		 //   std::cout << std::endl;
	//} 
	//int ioss;
	//std::cin >> ioss;

	//Remove duplicates去除重复项
	vector<NodePair> uniqueEdges(edges.size());
	size_t e = 0;
	NodePair last(numeric_limits<uint64_t>::max(), numeric_limits<uint64_t>::max());
	for (const NodePair& a : edges) {
		if (last.idA != a.idA || last.idB != a.idB) {
			last = a;
			uniqueEdges[e++] = NodePair(mapExternalNodeId(a.idA), a.idB);
		}
	}
	uniqueEdges.resize(e);
	for (NodePair& a : uniqueEdges) {
		a.idB = mapExternalNodeId(a.idB);
	}
	//LOG_PRINT("size() of uniqueEdges = " << uniqueEdges.size());
	//LOG_PRINT("[LOADING] Sorted edges");
	if (isUCG)
		LOG_PRINT("[LOADING] Number of nodes: " << nextNodeId);

	return GraphData(nextNodeId, move(uniqueEdges), move(revNodeRenaming), move(pro));
}*/

void GraphData::quick_sort(vector<NodePair>& edges,vector<float>& pro ,int l, int r)
{
	if (l < r)
	{
		//Swap(s[l], s[(l + r) / 2]); //将中间的这个数和第一个数交换 参见注1
		int i = l, j = r;
		auto x = edges[l];
		auto p = pro[l];
		while (i < j)
		{
			while (i < j && (edges[j].idA > x.idA ||(edges[j].idA==x.idA && edges[j].idB>x.idB))) // 从右向左找第一个小于x的数
				j--;
			if (i < j) {
				edges[i] = edges[j];
				pro[i++] = pro[j];
			}
			
			while (i < j && (edges[i].idA < x.idA || (edges[i].idA == x.idA && edges[i].idB < x.idB))) // 从左向右找第一个大于等于x的数
				i++;
			if (i < j) {
				edges[j] = edges[i];
				pro[j--] = pro[i];
			}
		}
		edges[i] = x;
		pro[i] = p;
		quick_sort(edges,pro, l, i - 1); // 递归调用 
		quick_sort(edges,pro, i + 1, r);
	}
}

GraphData GraphData::sampledata(const std::vector<NodePair>& GraphEdges, const std::vector<float>& GraphWeight) {//重载loadFromPath函数
	
	std::vector<NodePair> edges;
	//std::vector<float> pro;
	/*for (int i = 0; i < GraphEdges.size();i++) {
		int prob = GraphWeight[i] * BITYPE_WIDTH * sizeof(BITYPE) * 8;
		std::random_device rd;
		std::mt19937 mt(rd());
		std::uniform_int_distribution<> dis(1, BITYPE_WIDTH * sizeof(BITYPE) * 8);
		if (dis(mt) < prob) {
			continue;
		}
		//pro.push_back(GraphWeight[i]);
		edges.push_back(GraphEdges[i]);
	}*/

	std::unordered_map<uint64_t, uint64_t> nodeRenaming;
	std::unordered_map<uint64_t, uint64_t> revNodeRenaming;
	uint64_t nextNodeId = 0;

	auto mapExternalNodeId = [&nodeRenaming, &revNodeRenaming, &nextNodeId](uint64_t id) -> uint64_t {//Lambda表达式
		if (nodeRenaming.find(id) == nodeRenaming.end()) {
			nodeRenaming[id] = nextNodeId;//记录id对应的伪id
			revNodeRenaming[nextNodeId] = id;//记录伪id对应的真实id
			return nextNodeId++;
		}
		else {
			return nodeRenaming[id];//返回id对应的伪id
		}
	};


	//Remove duplicates去除重复项
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

/*GraphData GraphData::sampledata(const std::string& edgesFile, bool isUCG) {//重载loadFromPath函数
   // Count number of persons (excluding header line)
	std::string data_path = "../" + edgesFile;
	size_t numEdges = io::fileLines(data_path) - 1;

	const char *sep = "/";
	char path_str[100];
	strcpy(path_str, edgesFile.c_str());
	char *p;
	p = strtok(path_str, sep);
	p = strtok(NULL, sep);
	p = strtok(NULL, sep);
	string filename = p;
	std::string sample_data = "../test_queries/data/sample_" + filename;

	int fd = open(sample_data.c_str(), O_RDONLY);
	if (fd < 0) {
		FILE *fp;
		if ((fp = fopen(sample_data.c_str(), "w")) != NULL) {
			std::random_device rd;
			std::mt19937 g(rd());
			std::uniform_int_distribution<> dis(1, 100);
			for (int i = 0; i < numEdges; i++) {
				fprintf(fp, "%1.2f\n", float(dis(g)) / 100);
			}
		}
		fclose(fp);
		//新建一个sample_文件，每行写上概率。
	}
	io::MmapedFile file(data_path, O_RDONLY);

	Tokenizer tokenizer(file.mapping, file.size);
	tokenizer.skipLine(); // Skip header line 

	vector<NodePair> edges;
	edges.reserve(numEdges);
	vector<float> pro;
	pro.reserve(numEdges);
	io::MmapedFile sample_file(sample_data, O_RDONLY);
	Tokenizer sample_tokenizer(sample_file.mapping, sample_file.size);
	std::unordered_map<uint64_t, uint64_t> nodeRenaming;
	std::unordered_map<uint64_t, uint64_t> revNodeRenaming;
	uint64_t nextNodeId = 0;

	auto mapExternalNodeId = [&nodeRenaming, &revNodeRenaming, &nextNodeId](uint64_t id) -> uint64_t {//Lambda表达式
		if (nodeRenaming.find(id) == nodeRenaming.end()) {
			nodeRenaming[id] = nextNodeId;//记录id对应的伪id
			revNodeRenaming[nextNodeId] = id;//记录伪id对应的真实id
			return nextNodeId++;
		}
		else {
			return nodeRenaming[id];//返回id对应的伪id
		}
	};
	//vector<float>::iterator ite = pro.begin;
	while (!tokenizer.isFinished()&& !sample_tokenizer.isFinished()) {
		float tem_pro = atof((sample_tokenizer.readStr('\n')).c_str());
		int prob = tem_pro * BITYPE_WIDTH * sizeof(BITYPE) * 8;
		std::random_device rd;
		std::mt19937 mt(rd());
		std::uniform_int_distribution<> dis(1, BITYPE_WIDTH * sizeof(BITYPE) * 8);
		if (dis(mt) < prob) {
			sample_tokenizer.skipLine();
			tokenizer.skipLine();
			continue;
		}
		pro.push_back(tem_pro);

		assert(edges.size() < numEdges * 2);
		NodePair pair;
		pair.idA = tokenizer.readId('|');
		pair.idB = tokenizer.readId('\n');
		if (pair.idA == pair.idB) {
			continue; //No self-edges
		}

		//Add undirected（正向与反向）
		edges.push_back(pair);
		if (!isUCG)
			edges.push_back(NodePair(pair.idB, pair.idA));
	}

	
	//LOG_PRINT("size() of edges = " << edges.size());

	quick_sort(edges, pro, 0, pro.size()-1);

	/*for (int i = 0; i < edges.size(); i++) {
		std::cout << edges[i].idA << " " << edges[i].idB<<" | ";
		if (i % 20 == 0)
			std::cout << std::endl;
	} 


	vector<NodePair> uniqueEdges(edges.size());
	size_t e = 0;
	NodePair last(numeric_limits<uint64_t>::max(), numeric_limits<uint64_t>::max());
	for (const NodePair& a : edges) {
		if (last.idA != a.idA || last.idB != a.idB) {
			last = a;
			uniqueEdges[e++] = NodePair(mapExternalNodeId(a.idA), a.idB);
		}
	}
	uniqueEdges.resize(e);
	for (NodePair& a : uniqueEdges) {
		a.idB = mapExternalNodeId(a.idB);
	}
	//LOG_PRINT("size() of uniqueEdges = " << uniqueEdges.size());
	//LOG_PRINT("[LOADING] Sorted edges");
	if (isUCG)
		LOG_PRINT("[LOADING] Number of nodes: " << nextNodeId);
	std::vector<float> a; 
	return GraphData(nextNodeId, move(uniqueEdges), move(revNodeRenaming), a);
}*/

/*vector<NodePair> GraphData::loadComponent(const std::string& edgesFile) {//重载loadFromPath函数.
   // Count number of persons (excluding header line)
	size_t numEdges = io::fileLines(edgesFile) - 1;
	LOG_PRINT("[LOADING] Number of edges: " << numEdges);
	//cout << "test" << endl;

	// Load edges data from file
	LOG_PRINT("[LOADING] Loading edges from file: " << edgesFile);
	io::MmapedFile file(edgesFile, O_RDONLY);

	Tokenizer tokenizer(file.mapping, file.size);
	tokenizer.skipLine(); // Skip header line

	vector<NodePair> edges;
	edges.reserve(numEdges);

	std::unordered_map<uint64_t, uint64_t> nodeRenaming;
	std::unordered_map<uint64_t, uint64_t> revNodeRenaming;
	uint64_t nextNodeId = 0;

	auto mapExternalNodeId = [&nodeRenaming, &revNodeRenaming, &nextNodeId](uint64_t id) -> uint64_t {//Lambda表达式
		if (nodeRenaming.find(id) == nodeRenaming.end()) {
			nodeRenaming[id] = nextNodeId;//记录id对应的伪id
			revNodeRenaming[nextNodeId] = id;//记录伪id对应的真实id
			return nextNodeId++;
		}
		else {
			return nodeRenaming[id];//返回id对应的伪id.
		}
	};

	while (!tokenizer.isFinished()) {
		assert(edges.size() < numEdges * 2);
		NodePair pair;
		pair.idA = tokenizer.readId('|');
		pair.idB = tokenizer.readId('\n');
		if (pair.idA == pair.idB) {
			continue; //No self-edges
		}

		//Add undirected（正向与反向）
		edges.push_back(pair);
		edges.push_back(NodePair(pair.idB, pair.idA));
	}
	LOG_PRINT("size() of edges = " << edges.size()); 

	// Reading edge
	LOG_PRINT("[LOADING] Read edges");
	std::sort(edges.begin(), edges.end(), [](const NodePair& a, const NodePair& b) {//Lambda
		return a.idA < b.idA || (a.idA == b.idA && a.idB < b.idB);
	});

	//Remove duplicates去除重复项.
	vector<NodePair> uniqueEdges(edges.size());
	size_t e = 0;
	NodePair last(numeric_limits<uint64_t>::max(), numeric_limits<uint64_t>::max());
	for (const NodePair& a : edges) {
		if (last.idA != a.idA || last.idB != a.idB) {
			last = a;
			uniqueEdges[e++] = NodePair(mapExternalNodeId(a.idA), a.idB);
		}
	}
	uniqueEdges.resize(e);
	for (NodePair& a : uniqueEdges) {  
		a.idB = mapExternalNodeId(a.idB);
	}
	LOG_PRINT("size() of uniqueEdges = " << uniqueEdges.size());
	LOG_PRINT("[LOADING] Sorted edges");

	LOG_PRINT("[LOADING] Number of nodes: " << nextNodeId);

	return move(uniqueEdges);
}*/
