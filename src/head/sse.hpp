/**
Copyright (C) 2023/03/08 by Zhenfang Liu, Jianxiong Ye.
Code must not be used, distributed, without written consent by the authors.
*/
#pragma once

#include <immintrin.h>
#ifdef __AVX2__
#define AVX2 1
#endif

extern __m128i sseMasks[128];

#ifdef AVX2
extern __m256i sseMasks256[256];
#endif
