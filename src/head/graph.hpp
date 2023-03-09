/**
Copyright (C) 2023/03/08 by Zhenfang Liu, Jianxiong Ye.
Code must not be used, distributed, without written consent by the authors.
*/
#pragma once
#include"GraphFromFile.hpp"
#include "log.hpp"
#include "queue.hpp"
#include <iostream>
#include <cstdint>
#include <cstddef>
#include <cassert>
#include <utility>
#include <vector>
#include <string>
#include <limits>
#include <unordered_map>
#include <cstdlib>
#include <ctime>
#include "bitops.hpp"
//#include "idschedulers.hpp"
#include "sse.hpp"
#include <bitset>
#include <immintrin.h>
#include <emmintrin.h>
#include <smmintrin.h>
#include <random>
#include <algorithm>
//#include "base.hpp"

#define BITYPE  __m256i
#define BITYPE_WIDTH 4
#define SAMPLE_NUMS sizeof(BITYPE)*8*BITYPE_WIDTH
//Here SAMPLE_NUMS=1024


/*struct NodePair {
   uint64_t idA;		//uint64_t is long long
   uint64_t idB;

   NodePair() { }
   NodePair(uint64_t idA, uint64_t idB)
      : idA(idA), idB(idB)
   { }
};*/

struct GraphData {
   size_t numNodes;			//顶点数
   std::vector<NodePair> edges;//<0,100><0,104><1,233><1,999><2,3>...<100,0><104,0>...Here are virtual ID
   std::unordered_map<uint64_t,uint64_t> revNodeRenaming;//<0,v1><1,v2><2,v3>...record the real ID of virtual ID
   GraphData(const size_t numNodes, std::vector<NodePair> edges, std::unordered_map<uint64_t,uint64_t> revNodeRenaming)
      : numNodes(numNodes), edges(move(edges)), revNodeRenaming(std::move(revNodeRenaming)){
   }

   GraphData(GraphData& other) = delete;
   GraphData(GraphData&& other) = default;

   static GraphData loadFromPath(const std::vector<NodePair>& GraphEdges, const std::vector<float>& GraphWeight);
   static GraphData sampledata(const std::vector<NodePair>& GraphEdges, const std::vector<float>& GraphWeight);
   //static std::vector<NodePair> loadComponent(const std::string& edgesFile);
   static void quick_sort(std::vector<NodePair>& edges, std::vector<float>& pro, int l, int r);
};

template<class EntryType>
class SizedList {
public:
   typedef EntryType Size;
   typedef EntryType Entry;

private:
   Size count;

public:
   SizedList() {
   }

   /// Reads the size for the list at the pointer position
   inline const Size& size() const {
      return count;
   }

   inline Size& size() {
      return count;
   }

   /// Sets the size for the list at the pointer position
   void setSize(const Size count) {
      this->count = count;
   }

   const Entry* getPtr(const size_t index) const __attribute__ ((pure)) {
      const auto offsetPtr = reinterpret_cast<const uint8_t*>(this) + sizeof(Size) + sizeof(Entry)*index;
      return reinterpret_cast<const Entry*>(offsetPtr);
   }

   Entry* getPtr(const size_t index) __attribute__ ((pure)) {
      const auto offsetPtr = reinterpret_cast<uint8_t*>(this) + sizeof(Size) + sizeof(Entry)*index;
      return reinterpret_cast<Entry*>(offsetPtr);
   }

   std::pair<const Entry*,const Entry* const> bounds() const __attribute__ ((pure)) {
      const Entry* i=getPtr(0);
      return std::make_pair(i,i+count);
   }

   SizedList<Entry>* nextList(Size count) __attribute__ ((pure)) {
      auto offsetPtr = reinterpret_cast<uint8_t*>(this)+sizeof(Size)+sizeof(Entry)*count;
      return reinterpret_cast<SizedList<Entry>*>(offsetPtr);
   }
};



template<class EntryType>
class EdgesList {
public:
	typedef uint32_t Size;
	typedef EntryType Entry;

private:
	Size count;

public:
	EdgesList() {
	}

	/// Reads the size for the list at the pointer position
	inline const Size& size() const {
		return count;
	}

	inline Size& size() {
		return count;
	}

	/// Sets the size for the list at the pointer position
	void setSize(const Size count) {
		this->count = count;
	}

	const Entry* getPtr(const size_t index) const __attribute__((pure)) {

		const auto offsetPtr = reinterpret_cast<const uint8_t*>(this) + sizeof(Size) + sizeof(Entry)*index;
		return reinterpret_cast<const Entry*>(offsetPtr);
	}

	Entry* getPtr(const size_t index) __attribute__((pure)) {
		const auto offsetPtr = reinterpret_cast<uint8_t*>(this) + sizeof(Size) + sizeof(Entry)*index;
		return reinterpret_cast<Entry*>(offsetPtr);
	}

	std::pair<const Entry*, const Entry* const> bounds() const __attribute__((pure)) {
		const Entry* i = getPtr(0);
		return std::make_pair(i, i + count);
	}

	EdgesList<Entry>* nextList(Size count) __attribute__((pure)) {
		auto offsetPtr = reinterpret_cast<uint8_t*>(this) + sizeof(Size) + sizeof(Entry)*count;
		return reinterpret_cast<EdgesList<Entry>*>(offsetPtr);
	}
};


template<typename bit_t, uint64_t width>
struct EdgesBits {
	static const size_t TYPE_BITS_COUNT = sizeof(bit_t) * 8;

	bit_t data[width];

	EdgesBits() :data() {
		/*for (int i = 0; i < width; i++) {
			data[i] = BitBaseOp<bit_t>::zero();

		}*/

	}

	void getedges(float p) {
		int insert_num = 0;
		int insert_type = 0;
		if (p < 0.5) {
			insert_num = p * TYPE_BITS_COUNT*width;
			insert_type = 1;
		}
		else {
			insert_num = (1 - p) * TYPE_BITS_COUNT*width;
			insert_type = 0;
		}
		std::random_device rd;
		std::mt19937 mt(rd());
		std::uniform_int_distribution<> dis(1, TYPE_BITS_COUNT*width);
		//EdgesBits<bit_t,width> edge;
		if (insert_type == 1) {
			for (int i = 0; i < insert_num; i++) {
				int pos, field, field_bit;
				do {
					pos = dis(mt) - 1;
					field = pos / TYPE_BITS_COUNT;
					field_bit = pos - (field*TYPE_BITS_COUNT);
				} while (BitBaseOp<bit_t>::notZero(data[field] & BitBaseOp<bit_t>::getSetMask(field_bit)));
				data[field] |= BitBaseOp<bit_t>::getSetMask(field_bit);
			}

		}
		else {
			for (int i = 0; i < insert_num; i++) {
				int pos, field, field_bit;
				do {
					pos = dis(mt) - 1;
					field = pos / TYPE_BITS_COUNT;
					field_bit = pos - (field*TYPE_BITS_COUNT);
				} while (BitBaseOp<bit_t>::notZero(data[field] & BitBaseOp<bit_t>::getSetMask(field_bit)));
				data[field] |= BitBaseOp<bit_t>::getSetMask(field_bit);
			}

			for (int j = 0; j < width; j++) {
				data[j] = ~data[j];
			}
		}
	}
};

template<class IdType>
class Graph {
public:
   typedef IdType Id;
   typedef SizedList<Id>* Content;
   typedef uint32_t ComponentId;//The ID of the connecting component
   typedef uint64_t ComponentSize;

   const size_t numVertices;
   size_t numEdges;
   //std::vector<float> edges_prob;
   std::vector<ComponentId> personComponents;//Record the connected component ID for each vertex
   std::vector<ComponentSize> componentSizes;//Number of connected components
   std::vector<ComponentSize> componentEdgeCount;//
   ComponentSize maxComponentSize;

private:
   Content* table;

   std::unordered_map<uint64_t,uint64_t> revNodeRenaming;

public:
   uint8_t* data;
   EdgesBits<BITYPE,BITYPE_WIDTH>** edgesbit;
   //Constructor function, allocating numVertices size spaces for personComponents, each with 0
   Graph(size_t numVertices) : numVertices(numVertices), personComponents(numVertices), componentSizes(), componentEdgeCount(), maxComponentSize(), table(nullptr), data(nullptr), edgesbit(nullptr){
      // XXX: Maybe posix_memalign helps?
      table = new Content[numVertices]();
   }

   Graph(Graph& other) = delete;

   Graph(Graph&& other) : numVertices(other.numVertices), numEdges(other.numEdges), personComponents(other.personComponents), componentSizes(other.componentSizes), componentEdgeCount(other.componentEdgeCount), maxComponentSize(other.maxComponentSize), table(other.table), revNodeRenaming(std::move(other.revNodeRenaming)), data(other.data),edgesbit(other.edgesbit) {
      other.table=nullptr;
      other.data=nullptr;
	  other.edgesbit = nullptr;
   }

   ~Graph() {
      if(table) {
         delete[] table;
         table = nullptr;
      }

      if(data) {
         delete[] data;
         data = nullptr;
      }

	  if (edgesbit) {
		  delete[] edgesbit;
		  edgesbit = nullptr;
	  }

   }

   IdType mapInternalNodeId(IdType id) const {
      const auto iter = revNodeRenaming.find(id);
      if(iter!=revNodeRenaming.cend()) {
         return iter->second;
      } else {
         throw -1;
      }
   }

   /// Inserts the data for the specified id into the index
   void insert(Id id, Content content) {
      assert(table!=nullptr);
      assert(id<numVertices);
      table[id] = std::move(content);
   }

   /// Retrieves the data for the specified id
   __attribute__ ((pure)) Content retrieve(Id id) const {
      assert(id<numVertices);
      assert(table[id]!=nullptr);
      return table[id];
   }

   /*__attribute__((pure)) BITYPE findedges(Id id, int friendId, uint32_t width) const{
	   return edgesbit[id][friendId].data[width];
   }*/

   inline IdType maxKey() const {
      return numVertices-1;
   }

   inline IdType size() const {
      return numVertices;
   }

   static Graph sampleGraph(const std::vector<NodePair>& GraphEdges,const std::vector<float>& GraphWeight) {//Generate a persongraph from graphdata, and the table contains neighbor nodes
	   GraphData graphData = GraphData::sampledata(GraphEdges, GraphWeight);
	   IdType numPersons = graphData.numNodes;//Number of vertices
	   std::vector<NodePair>& edges = graphData.edges;//Number of edges
	   Graph personGraph(numPersons);
	   //personGraph.edges_prob = std::move(graphData.prob);
	   personGraph.revNodeRenaming = std::move(graphData.revNodeRenaming);

	   std::vector<NodePair> double_edges;//Bidirectional edges are used to divide the connected components
	   double_edges.reserve(edges.size() * 2);
	   for (NodePair& a : edges) {
		   double_edges.push_back(a);
		   double_edges.push_back(NodePair(a.idB, a.idA));
	   }
	   std::sort(double_edges.begin(), double_edges.end(), [](const NodePair& a, const NodePair& b) {//Lambda
		   return a.idA < b.idA || (a.idA == b.idA && a.idB < b.idB);
	   });
	   std::vector<NodePair> uniqueEdges(double_edges.size());
	   size_t e = 0;
	   std::vector<size_t> index;
	   NodePair last(std::numeric_limits<uint64_t>::max(), std::numeric_limits<uint64_t>::max());
	   for (const NodePair& a : double_edges) {
		   if (last.idA != a.idA) {
			   index.push_back(e);
			   last = a;
			   uniqueEdges[e++] = NodePair(a.idA, a.idB);
		   }
		   else if (last.idB != a.idB) {
			   last = a;
			   uniqueEdges[e++] = NodePair(a.idA, a.idB);
		   }
	   }
	   uniqueEdges.resize(e);
	   index.push_back(e);
	   const size_t dataSize = (numPersons + edges.size()) * sizeof(IdType);
	   //std::cout << "start sample one graph" << std::endl;
	   uint8_t* data = new uint8_t[dataSize]();
	   {//Add neighbors for each vertexof the personGraph
		   SizedList<IdType>* neighbours = reinterpret_cast<SizedList<IdType>*>(data);
		   size_t ix = 0;

		   std::random_device rd;
		   std::mt19937 mt(rd());
		   std::uniform_int_distribution<> dis(1, 100);
		   for (IdType person = 0; person < numPersons; person++) {
			   assert(reinterpret_cast<uint8_t*>(neighbours) < data + dataSize);
			   IdType* insertPtr = neighbours->getPtr(0);
			   IdType count = 0;
			   while (ix < edges.size() && edges[ix].idA == person) {
				   int prob = GraphWeight[ix] * 100;
				   int pr = dis(mt);
				   if (pr > prob) {
					   ix++;
					   continue;
				   }
				   *insertPtr = edges[ix].idB;
				   insertPtr++;
				   count++;
				   ix++;
			   }

			   neighbours->setSize(count);
			   personGraph.insert(person, neighbours);
			   neighbours = neighbours->nextList(count);
		   }
	   }

	   personGraph.data = data;
	   personGraph.numEdges = edges.size();

	   //LOG_PRINT("[LOADING] Created person graph of size: " << dataSize / 1024 << " kb");
	   personGraph.analyzeGraph(uniqueEdges,index);

	   return std::move(personGraph);
   }

   static Graph loadUCGFromPath(const std::vector<NodePair>& GraphEdges, const std::vector<float>& GraphWeight) {//Generate a persongraph from graphdata, and the table contains neighbor nodes
	   //LOG_PRINT("LoadUCGFromPath");
	   GraphData graphData = GraphData::loadFromPath(GraphEdges, GraphWeight);
	   //std::vector<NodePair> tempedges = std::move(GraphData::loadComponent(edgesFile))
	   IdType numPersons = graphData.numNodes;//number of vertices
	   std::vector<NodePair>& edges = graphData.edges;//number of edges
	   //LOG_PRINT("Nodes: " << numPersons);
	  // LOG_PRINT("Edges: " << edges.size());

	   std::vector<NodePair> double_edges;
	   double_edges.reserve(edges.size() * 2);
	   for (NodePair& a : edges) {
		   double_edges.push_back(a);
		   double_edges.push_back(NodePair(a.idB, a.idA));
	   }
	   std::sort(double_edges.begin(), double_edges.end(), [](const NodePair& a, const NodePair& b) {//Lambda
		   return a.idA < b.idA || (a.idA == b.idA && a.idB < b.idB);
	   });
	   std::vector<NodePair> uniqueEdges(double_edges.size());
	   size_t e = 0;
	   std::vector<size_t> index;
	   NodePair last(std::numeric_limits<uint64_t>::max(), std::numeric_limits<uint64_t>::max());
	   for (const NodePair& a : double_edges) {
		   if (last.idA != a.idA) {
			   index.push_back(e);
			   last = a;
			   uniqueEdges[e++] = NodePair(a.idA, a.idB);
		   }
		   else if (last.idB != a.idB) {
			   last = a;
			   uniqueEdges[e++] = NodePair(a.idA, a.idB);
		   }
	   }
	   uniqueEdges.resize(e);
	   index.push_back(e);

	  // for (int i = 0; i < 500; i++) {
		   //std::cout << uniqueEdges[i].idA << "--" << uniqueEdges[i].idB <<" "<< index[uniqueEdges[i].idA] << std::endl;
	  // }
	   //LOG_PRINT("Graph personGraph(numPersons)");
	   Graph personGraph(numPersons);
	   //LOG_PRINT("after Graph personGraph(numPersons)");
	   personGraph.revNodeRenaming = std::move(graphData.revNodeRenaming);
	   const size_t dataSize = (numPersons + edges.size()) * sizeof(IdType);
	   uint8_t* data = new uint8_t[dataSize]();

	   {//为personGraph的每个点添加邻居
		   SizedList<IdType>* neighbours = reinterpret_cast<SizedList<IdType>*>(data);
		   size_t ix = 0;
		   for (IdType person = 0; person < numPersons; person++) {
			   assert(reinterpret_cast<uint8_t*>(neighbours) < data + dataSize);
			   IdType* insertPtr = neighbours->getPtr(0);
			   IdType count = 0;
			   while (ix < edges.size() && edges[ix].idA == person) {
				   *insertPtr = edges[ix].idB;
				   insertPtr++;
				   count++;
				   ix++;
			   }
			   neighbours->setSize(count);
			   personGraph.insert(person, neighbours);
			   neighbours = neighbours->nextList(count);
		   }
	   }

	   personGraph.data = data;
	   personGraph.numEdges = edges.size();
	   EdgesBits<BITYPE, BITYPE_WIDTH>** temp_edge;
	   temp_edge = new EdgesBits<BITYPE, BITYPE_WIDTH>*[numPersons]();
	   for (int a = 0; a < numPersons; a++) {
		   const auto& curFriends = *personGraph.retrieve(a);
		   int size = curFriends.size();
		   //std::cout << personGraph.mapInternalNodeId(a) << "|" << size << " ";
		   const auto ret = posix_memalign(reinterpret_cast<void**>(&(temp_edge[a])), 64, sizeof(EdgesBits<BITYPE, BITYPE_WIDTH>)*size);//Data alignment (allocate first memory address, align boundary, specify allocation byte size)

		   if (unlikely(ret != 0)) {

			   std::cout << "unlikely" << std::endl;
			   throw - 1;
		   }

		   new(temp_edge[a]) EdgesBits<BITYPE, BITYPE_WIDTH>[size]();
	   }
	   personGraph.edgesbit = std::move(temp_edge);
	   {
		   size_t iy = 0;
		   size_t edgeid = 0;
		   for (IdType person = 0; person < numPersons; person++) {
			   IdType count = 0;
			   while (iy < edges.size() && edges[iy].idA == person) {
				   personGraph.edgesbit[person][count].getedges(GraphWeight[edgeid]);
				   int numof1 = 0;
				   for (int i = 0; i < BITYPE_WIDTH; i++){
					   numof1 += BitBaseOp<BITYPE>::popCount(personGraph.edgesbit[person][count].data[i]);
					}
				   //std::cout << personGraph.mapInternalNodeId(person)<<" " << personGraph.mapInternalNodeId(edges[iy].idB)<<" "<< numof1 << std::endl;

				   count++;
				   iy++;
				   edgeid++;
			   }
		   }
	   }


	  /*const auto& curFriendsedges = personGraph.retrieveedges(0);
	   auto friendsedgesBounds = curFriendsedges->bounds();
	   while (friendsedgesBounds.first != friendsedgesBounds.second) {
		   std::cout << (friendsedgesBounds.first) << "address";
		   std::cout << BitBaseOp<BITYPE>::popCount(friendsedgesBounds.first->data[0]) << std::endl;

		   friendsedgesBounds.first++;
	  }*/

		//LOG_PRINT("[LOADING] Created person graph of size: " << dataSize / 1024 << " kb");

		//LOG_PRINT("start analyze");
		personGraph.analyzeGraph(uniqueEdges,index);
		return std::move(personGraph);
   }

private:
	void analyzeGraph(std::vector<NodePair>& edges, std::vector<size_t>& index) {//Divide the connectivity component
		const auto graphSize = size();

		componentSizes.push_back(std::numeric_limits<ComponentSize>::max()); // Component 0 is invalid
		componentEdgeCount.push_back(std::numeric_limits<ComponentSize>::max()); // Component 0 is invalid

		// Identify connected components by running bfs
		awfy::FixedSizeQueue<IdType> toVisit = awfy::FixedSizeQueue<IdType>(graphSize);

		size_t trivialComponents = 0;
		ComponentId componentId = 1;//Connecting component IDs, 1, 2, 3...
		for (IdType node = 0; node < graphSize; node++) {
			if (personComponents[node] != 0) { continue; }//personComponents initializes with a default value of 0

			ComponentSize componentSize = 1;
			personComponents[node] = componentId;
			toVisit.push_back_pos() = node;

			uint64_t componentNeighborCount = 0;

			// Do BFS for one connected component
			do {
				const IdType curNode = toVisit.front();
				toVisit.pop_front();


				/*for (const NodePair& a : edges) {
					if (a.idA > curNode)
						break;
					else if (a.idA < curNode)
						continue;
					else {
						if (personComponents[a.idB] != 0)
							continue;
						personComponents[a.idB] = componentId;
						componentSize++;
						toVisit.push_back_pos() = a.idB;
						componentNeighborCount++;
					}
				}*/

				for (int i = index[curNode]; i < index[curNode + 1]; i++) {
					if (personComponents[edges[i].idB] != 0)
						continue;
					personComponents[edges[i].idB] = componentId;
					componentSize++;
					toVisit.push_back_pos() = edges[i].idB;
					componentNeighborCount++;
				}

				/*const auto curNeighbours = retrieve(curNode);
				componentNeighborCount += curNeighbours->size();//Record the number of edges for each connected component (one edge is recorded twice - front and back)
				auto neighbourBounds = curNeighbours->bounds();//<i,j> i is the vertex id(0,1,2...) The pointer, j is the pointer to the number of vertex neighbors +i
				while (neighbourBounds.first != neighbourBounds.second) {
					const IdType curNeighbour = *neighbourBounds.first;
					++neighbourBounds.first;
					if (personComponents[curNeighbour] != 0) { continue; }//Duplicates are excluded, such as 0 and 999 neighbors, and 0 is not repeatedly inserted into toVisit
					personComponents[curNeighbour] = componentId;
					componentSize++;
					toVisit.push_back_pos() = curNeighbour;
				}*/
			} while (!toVisit.empty());

			componentEdgeCount.push_back(componentNeighborCount); //Record the number of edges corresponding to the connected component
			componentSizes.push_back(componentSize);//Record the number of vertices corresponding to the connected component
			if (componentSize > maxComponentSize) {//maxComponentSize records the number of vertices of the largest connected component (the more vertices, the larger the connected component)
				maxComponentSize = componentSize;
			}
			if (componentSize < 5) {
				trivialComponents++;
			}
			else {
				//std::cout << "# C " << componentSize << std::endl;
			}

			componentId++;

			// Check for overflow
			assert(componentId > 0);
		}

		//LOG_PRINT("[Query4] Max component size " << maxComponentSize);
		//std::cout << "# Found number components  " << componentId - 1 << " (" << trivialComponents << " are of size < 5)." << std::endl;
	}

   void analyzeGraph() {//划分连通分量
      const auto graphSize = size();

      componentSizes.push_back(std::numeric_limits<ComponentSize>::max()); // Component 0 is invalid
      componentEdgeCount.push_back(std::numeric_limits<ComponentSize>::max()); // Component 0 is invalid

      // Identify connected components by running bfs
      awfy::FixedSizeQueue<IdType> toVisit = awfy::FixedSizeQueue<IdType>(graphSize);

      size_t trivialComponents=0;
      ComponentId componentId=1;//Connecting component IDs, 1, 2, 3...
      for(IdType node=0; node<graphSize; node++) {
         if(personComponents[node]!=0) { continue; }//personComponents initializes with a default value of 0

         ComponentSize componentSize=1;
         personComponents[node]=componentId;
         toVisit.push_back_pos()=node;

         uint64_t componentNeighborCount=0;

         // Do BFS for one connected component
         do {
            const IdType curNode = toVisit.front();
            toVisit.pop_front();

            const auto curNeighbours=retrieve(curNode);
            componentNeighborCount += curNeighbours->size();//Record the number of edges for each connected component (one edge is recorded twice - front and back)

            auto neighbourBounds = curNeighbours->bounds();//<i,j>  i is the vertex id(0,1,2...) The pointer, j is the pointer to the number of vertex neighbors +i
            while(neighbourBounds.first != neighbourBounds.second) {
               const IdType curNeighbour=*neighbourBounds.first;
               ++neighbourBounds.first;
               if (personComponents[curNeighbour]!=0) { continue; }//Duplicates are excluded, such as 0 and 999 neighbors, and 0 is not repeatedly inserted into toVisit
               personComponents[curNeighbour]=componentId;
               componentSize++;
               toVisit.push_back_pos() = curNeighbour;
            }
         } while(!toVisit.empty());

         componentEdgeCount.push_back(componentNeighborCount); //Record the number of edges corresponding to the connected component
         componentSizes.push_back(componentSize);//Record the number of vertices corresponding to the connected component
         if(componentSize>maxComponentSize) {//maxComponentSize records the number of vertices of the largest connected component (the more vertices, the larger the connected component)
            maxComponentSize = componentSize;
         }
         if(componentSize<5) {
            trivialComponents++;
         }
         componentId++;

         // Check for overflow
         assert(componentId>0);
      }

	  /*for (int i = 0; i < graphSize; i++) {
		  std::cout << "node " << i << " 's id= " << personComponents[i] << std::endl;
	  }*/

	  //std::cout << "maxComponentSize = " << maxComponentSize << ", num of components = " << componentId - 1 << ", maxComponentSize = " << maxComponentSize << std::endl;
	  //for (int i = 1; i < componentId ; i++) {
		 // std::cout << "componentid " << i << " has " << componentSizes[i] << " nodes" << std::endl;
	  //}

      //LOG_PRINT("[Query4] Max component size "<< maxComponentSize);
      //std::cout<<"# Found number components  "<< componentId-1<<" ("<<trivialComponents<<" are of size < 5)."<<std::endl;
   }
};
