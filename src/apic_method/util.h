#pragma once

#include <algorithm>
#include <cmath>
#include <vector>


#if defined(USE_TBB)
#include <tbb/tbb.h>
#elif !defined(SINGLE_THREADED)
#include <execution>
#include <thread>
#endif

#include "math_defs.h"

using std::max;
using std::min;

/// Parallel for loop with TBB or C++17 fallback
template <typename Index, typename Callable>
inline void parallel_for(Index start, Index end, Callable c) {
#if defined(USE_TBB)
  tbb::parallel_for(start, end, [&](Index i) { c(i); });
#elif !defined(SINGLE_THREADED)
  static const Index n_threads = std::thread::hardware_concurrency();
  static std::vector<Index> worker_pool(n_threads);

  const Index work_per_thread = (end - start + n_threads - 1) / n_threads;

  for (Index thread_id = 0; thread_id < n_threads; ++thread_id) {
    for (Index offset = 0; offset < work_per_thread; ++offset) {
      Index i = start + thread_id * work_per_thread + offset;
      if (i >= end)
        continue;
      c(i);
    }
  }
#else
  for (Index i = start; i < end; ++i) {
    c(i);
  }
#endif
}

/// Square a value
template <typename T>
inline T sqr(const T& x) {
  return x * x;
}

/// Clamp a value within a specified range
template <typename T>
inline T clamp(T value, T lower, T upper) {
  if (value < lower) return lower;
  if (value > upper) return upper;
  return value;
}

/// Compute cell index and barycentric coordinate
template <typename T>
inline void get_barycentric(T x, int& base, T& fraction, int lower_bound, int upper_bound) {
  T floored = std::floor(x);
  base = static_cast<int>(floored);

  if (base < lower_bound) {
    base = lower_bound;
    fraction = T(0);
  } else if (base > upper_bound - 2) {
    base = upper_bound - 2;
    fraction = T(1);
  } else {
    fraction = x - floored;
  }
}

/// Zero out a std::vector
template <typename T>
void zero(std::vector<T>& vec) {
  for (int i = static_cast<int>(vec.size()) - 1; i >= 0; --i) {
    vec[i] = T(0);
  }
}

/// Insert an element at a given index by shifting right
template <typename T>
void insert(std::vector<T>& arr, unsigned int index, T val) {
  arr.push_back(arr.back());
  for (int i = static_cast<int>(arr.size()) - 1; i > static_cast<int>(index); --i) {
    arr[i] = arr[i - 1];
  }
  arr[index] = val;
}