/**
Copyright (C) 2023/03/08 by Zhenfang Liu, Jianxiong Ye.
Code must not be used, distributed, without written consent by the authors.
*/
#pragma once

#define likely(x)       __builtin_expect(!!(x), 1)
#define unlikely(x)     __builtin_expect(!!(x), 0)//Reduce jumps because ifelse is biased

#ifndef GCC_VERSION
#define __builtin_assume_aligned(a, b) a
#endif
