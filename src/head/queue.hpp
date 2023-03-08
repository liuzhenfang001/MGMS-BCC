/**
Copyright (C) 2023/03/08 by Zhenfang Liu, Jianxiong Ye.
Code must not be used, distributed, without written consent by the authors.
*/
#pragma once

#include <cstddef>
#include <cassert>
#include <utility>
#include "macros.hpp"
#include <iostream>
namespace awfy {

   template<class T>
   class FixedSizeQueue {
   private:
	   T* elems;
	   T* startPtr;
	   T* endPtr;
	   size_t size_;

   public:
	   FixedSizeQueue(size_t size) : elems(new T[size]), startPtr(elems), endPtr(elems), size_(size)
	   {
	   }

	   ~FixedSizeQueue() {
		   if (elems != nullptr) {
			   delete[] elems;
		   }
	   }

	   inline bool empty() const {
		   return startPtr == endPtr;
	   }

	   inline size_t size() const {
		   return endPtr - startPtr;
	   }

	   inline bool isfull() const {
		   return endPtr < elems + size_ ? false : true;
	   }

	   inline const T& front() const {
		   assert(!empty());
		   return *startPtr;
	   }

	   inline void pop_front() {
		   assert(!empty());
		   startPtr++;
	   }

	   inline T& push_back_pos() {
		   assert(endPtr < elems + size_);
		   return *endPtr++;
	   }

	   void reset(size_t newSize) {
		   if (newSize > size_) {
			   delete[] elems;
			   elems = new T[newSize];
			   size_ = newSize;
		   }
		   startPtr = elems;
		   endPtr = elems;
	   }

	   inline std::pair<T*, T*> bounds() {
		   return make_pair(startPtr, endPtr);
	   }
   };
   /*class FixedSizeQueue {
   private:
      T* elems;
      T* startPtr;
      T* endPtr;
      size_t size_;

   public:
      FixedSizeQueue(size_t size) : elems(new T[size]),startPtr(elems),endPtr(elems),size_(size)
      {  }

      ~FixedSizeQueue() {
         if(elems!=nullptr) {
            delete[] elems;
         }
      }

      inline bool empty() const {
         return startPtr==endPtr;
      }

      inline size_t size() const {
		  return endPtr > startPtr ? endPtr - startPtr : endPtr + size_ - startPtr;
      }

      inline const T& front() const {
         assert(!empty());
         return *startPtr;
      }

      inline void pop_front() {
         assert(!empty());
		 //startPtr = elems + (startPtr - elems + 1) % size_;
		 if (startPtr == elems + size_) {
			 startPtr = elems;
		 }
		 else {
			 startPtr++;
		 }
      }

      inline T& push_back_pos() {
		 assert(endPtr != startPtr - 1);
		 if (endPtr < elems+size_) {
			 return *endPtr++;
		 }
		 else {
			 endPtr = elems;
			 return *(elems + size_);
		 }
      }

      void reset(size_t newSize) {
         if(newSize>size_) {
            delete[] elems;
            elems = new T[newSize];
            size_ = newSize;
         }
         startPtr=elems;
         endPtr=elems;
      }

      inline std::pair<T*,T*> bounds() {
         return make_pair(startPtr,endPtr);
      }
   };*/

   /*template<class T>
   class FixedSizeQueue1 {
   private:
	   T* elems;
	   T* startPtr;
	   T* endPtr;
	   size_t size_;

   public:
	   FixedSizeQueue1(size_t size)
	   {
		   const auto ret = posix_memalign(reinterpret_cast<void**>(&(elems)), 64, sizeof(T)*size);//
		   if (unlikely(ret != 0)) {
			   throw - 1;
		   }
		   new(elems) T[size]();
		   startPtr = elems;
		   endPtr = elems;
		   size_ = size;
	   }

	   ~FixedSizeQueue1() {
		   if (elems != nullptr) {
			   delete[] elems;
		   }
	   }

	   inline bool empty() const {
		   return startPtr == endPtr;
	   }

	   inline size_t size() const {
		   return endPtr - startPtr;
	   }

	   inline const T& front() const {
		   assert(!empty());
		   return *startPtr;
	   }

	   inline void pop_front() {
		   assert(!empty());
		   startPtr++;
	   }

	   inline T& push_back_pos() {
		   assert(endPtr < elems + size_);
		   return *endPtr++;
	   }

	   void reset(size_t newSize) {
		   if (newSize > size_) {
			   delete[] elems;
			   elems = new T[newSize];
			   size_ = newSize;
		   }
		   startPtr = elems;
		   endPtr = elems;
	   }

	   inline std::pair<T*, T*> bounds() {
		   return make_pair(startPtr, endPtr);
	   }
   };*/

   /*
   template<class T>
   class QuadrupleQueue {
   private:
	   T* elems;
	   T* startPtr;
	   T* endPtr;
	   size_t size_;

   public:
	   QuadrupleQueue(size_t size) : elems(new T[size]), startPtr(elems), endPtr(elems), size_(size)
	   {  }

	   ~QuadrupleQueue() {
		   if (elems != nullptr) {
			   delete[] elems;
		   }
	   }

	   inline bool empty() const {
		   return startPtr == endPtr;
	   }

	   inline bool isfull() const {
		   return   endPtr - startPtr == -1;
	   }

	   inline size_t size() const {
		   return ((endPtr - startPtr) + size_)% size_;
	   }

	   inline const T& front() const {
		   assert(!empty());
		   return *startPtr;
	   }

	   inline void pop_front() {
		   assert(!empty());
		   if (startPtr < elems + size_ - 1)
			   startPtr++;
		   else
			   startPtr = elems;
	   }

	   inline T& push_back_pos() {
		   //assert(endPtr < elems + size_);
		   if (isfull()){
			   FATAL_ERROR("queue has full when push_back_pos");
			   exit(-1);
		   }
		   else {
			   return *endPtr++;
		   }

	   }

	   void reset(size_t newSize) {
		   if (newSize > size_) {
			   delete[] elems;
			   elems = new T[newSize];
			   size_ = newSize;
		   }
		   startPtr = elems;
		   endPtr = elems;
	   }

	   inline std::pair<T*, T*> bounds() {
		   return make_pair(startPtr, endPtr);
	   }
   };

   template<class T>
   class QuadrupleQueue1 {
   private:
	   T* elems;
	   T* startPtr;
	   T* endPtr;
	   size_t size_;
	   QuadrupleQueue1<T>* next;
   public:
	   QuadrupleQueue1(size_t size) : elems(new T[size]), startPtr(elems), endPtr(elems), size_(size), next(NULL)
	   {  }

	   ~QuadrupleQueue1() {
		   if (elems != nullptr) {
			   delete[] elems;
			   next = NULL;
		   }
	   }

	   inline bool empty() const {
		   return (startPtr == endPtr) && (next != NULL ? next->empty() : true);
	   }

	   inline size_t size() const {
			   return (endPtr - startPtr) + next != NULL ? next->size() : 0;
	   }

	   inline const T& front() const {
		   assert(!empty());
		   if (startPtr == endPtr)
			   return next->front();
		   return *startPtr;
	   }

	   inline QuadrupleQueue<T> pop_front() {
		   assert(!empty());
		   if (startPtr == endPtr) {
				next->pop_front();
				return *next;
		   }
		   else {
				startPtr++;
				return this;
		   }

	   }

	   inline T& push_back_pos() {
		   assert(endPtr < elems + size_);
		   return *endPtr++;
	   }

	   void reset(size_t newSize) {
		   if (newSize > size_) {
			   delete[] elems;
			   elems = new T[newSize];
			   size_ = newSize;
		   }
		   startPtr = elems;
		   endPtr = elems;
	   }

	   inline std::pair<T*, T*> bounds() {
		   return make_pair(startPtr, endPtr);
	   }
   };
*/
}
