/**
Copyright (C) 2023/03/08 by Zhenfang Liu, Jianxiong Ye.
Code must not be used, distributed, without written consent by the authors.
*/
#include "head/GraphFromFile.hpp"
#include <iostream>
//#define unDIRECTED
void GraphFromFile::init(const std::string& edgesFile) {
	std::string data_path = "../" + edgesFile;
	size_t numEdges = io::fileLines(data_path) - 1;
	//std::cout << "Edges Lines: " << numEdges << std::endl;
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
			for (int i = 0; i < 2*numEdges; i++) {
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

	std::vector<NodePair> temp_edges;
	temp_edges.reserve(2*numEdges);

	while (!tokenizer.isFinished()) {
		assert(temp_edges.size() < numEdges * 3);
		NodePair pair;
		pair.idA = tokenizer.readId('|');
		pair.idB = tokenizer.readId('\n');
		if (pair.idA == pair.idB) {
			continue; //No self-edges
		}
		temp_edges.push_back(pair);
#ifdef unDIRECTED
	temp_edges.push_back(NodePair(pair.idB,pair.idA));
#endif // unDIRECTED


	}

	//quick_sort(temp_edges, pro, 0, numEdges - 1);
	std::sort(temp_edges.begin(), temp_edges.end(), [](const NodePair& a, const NodePair& b) {//Lambda
	   return a.idA<b.idA||(a.idA==b.idA && a.idB<b.idB);
	});
	//std::cout << "finish sort" << std::endl;

	//Remove duplicates去除重复项
	std::vector<NodePair> uniqueEdges(temp_edges.size());
	size_t e = 0;
	NodePair last(std::numeric_limits<uint64_t>::max(), std::numeric_limits<uint64_t>::max());
	for (const NodePair& a : temp_edges) {
		if (last.idA != a.idA || last.idB != a.idB) {
			last = a;
			uniqueEdges[e++] = NodePair(a.idA, a.idB);
		}
	}
	uniqueEdges.resize(e);
	//std::cout << "Edges Lines: " << e << std::endl;
	edges = std::move(uniqueEdges);
	weight = std::move(pro);
}

void GraphFromFile::getSource(const int maxBfs) {
}
