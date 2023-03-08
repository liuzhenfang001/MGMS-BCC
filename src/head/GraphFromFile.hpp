/**
Copyright (C) 2023/03/08 by Zhenfang Liu, Jianxiong Ye.
Code must not be used, distributed, without written consent by the authors.
*/
#pragma once
#include <stdlib.h>
#include <algorithm>
#include <string.h>
#include <random>
#include <vector>
#include "io.hpp"
#include "tokenizer.hpp"
struct NodePair {
	uint64_t idA;		//uint64_t is long long
	uint64_t idB;

	NodePair() { }
	NodePair(uint64_t idA, uint64_t idB)
		: idA(idA), idB(idB)
	{ }
};

class GraphFromFile {
public:
	std::vector<float> weight;
	std::vector<NodePair> edges;

public:
	GraphFromFile() :weight(),edges(){
	}
	void getSource(const int maxBfs);
	void init(const std::string& edgesFile);

	/*static void init(const std::string& edgesFile) {
		std::string data_path = "../" + edgesFile;
		size_t numEdges = io::fileLines(data_path) - 1;

		const char *sep = "/";
		char path_str[100];
		strcpy(path_str, edgesFile.c_str());
		char *p;
		p = strtok(path_str, sep);
		p = strtok(NULL, sep);
		p = strtok(NULL, sep);
		std::string filename = p;
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

		std::vector<float> pro;
		pro.reserve(numEdges);
		io::MmapedFile sample_file(sample_data, O_RDONLY);
		Tokenizer sample_tokenizer(sample_file.mapping, sample_file.size);
		while (!sample_tokenizer.isFinished()) {
			assert(pro.size() < numEdges * 2);
			float tem_pro = atof((sample_tokenizer.readStr('\n')).c_str());
			pro.push_back(tem_pro);
		}

		io::MmapedFile file(data_path, O_RDONLY);

		Tokenizer tokenizer(file.mapping, file.size);
		tokenizer.skipLine(); // Skip header line

		std::vector<NodePair> edges;
		edges.reserve(numEdges);

		while (!tokenizer.isFinished()) {
			assert(edges.size() < numEdges * 2);
			NodePair pair;
			pair.idA = tokenizer.readId('|');
			pair.idB = tokenizer.readId('\n');
			if (pair.idA == pair.idB) {
				continue; //No self-edges
			}
			edges.push_back(pair);
		}

		quick_sort(edges, pro, 0, numEdges - 1);

		//Remove duplicates去除重复项
		std::vector<NodePair> uniqueEdges(edges.size());
		size_t e = 0;
		NodePair last(std::numeric_limits<uint64_t>::max(), std::numeric_limits<uint64_t>::max());
		for (const NodePair& a : edges) {
			if (last.idA != a.idA || last.idB != a.idB) {
				last = a;
				uniqueEdges[e++] = NodePair(a.idA, a.idB);
			}
		}
		uniqueEdges.resize(e);

		GraphFromFile::edges = std::move(uniqueEdges);
		GraphFromFile::weight = std::move(pro);
	}*/

	static void quick_sort(std::vector<NodePair>& edges, std::vector<float>& pro, int l, int r){
		if (l < r){
			int i = l, j = r;
			auto x = edges[l];
			auto p = pro[l];
			while (i < j) {
				while (i < j && (edges[j].idA > x.idA || (edges[j].idA == x.idA && edges[j].idB > x.idB)))
					j--;
				if (i < j) {
					edges[i] = edges[j];
					pro[i++] = pro[j];
				}
				while (i < j && (edges[i].idA < x.idA || (edges[i].idA == x.idA && edges[i].idB < x.idB)))
					i++;
				if (i < j) {
					edges[j] = edges[i];
					pro[j--] = pro[i];
				}
			}
			edges[i] = x;
			pro[i] = p;
			quick_sort(edges, pro, l, i - 1);
			quick_sort(edges, pro, i + 1, r);
		}
	}
};
