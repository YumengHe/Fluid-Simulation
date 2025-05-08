#pragma once

#include <algorithm>
#include <cassert>
#include <climits>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "math_defs.h"


// This header provides:
// - Array1<T>: A dynamic, resizable 1D array class (for POD types)
// - WrapArray1<T>: A lightweight wrapper for externally managed 1D buffers
//
// Both follow STL-like semantics. WrapArray1 does not allocate memory or manage ownership.

struct Array1False {};  // Used to identify non-integral types in templates

template <typename T>
struct Array1IsIntegral {
  using type = Array1False;
};

// Array1<T>: Dynamic array for fundamental data types
template <typename T>
struct Array1 {
  using iterator = T*;
  using const_iterator = const T*;
  using size_type = unsigned long;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  unsigned long n = 0;
  unsigned long max_n = 0;
  T* data = nullptr;

  Array1() = default;

  Array1(const Array1<T>& x) {
    data = static_cast<T*>(std::malloc(x.n * sizeof(T)));
#ifndef NO_THROWS
    if (!data) throw std::bad_alloc();
#endif
    n = max_n = x.n;
    std::memcpy(data, x.data, n * sizeof(T));
  }

  Array1<T>& operator=(const Array1<T>& x) {
    if (max_n < x.n) {
      T* tmp = static_cast<T*>(std::malloc(x.n * sizeof(T)));
#ifndef NO_THROWS
      if (!tmp) throw std::bad_alloc();
#endif
      std::free(data);
      data = tmp;
      max_n = x.n;
    }
    n = x.n;
    std::memcpy(data, x.data, n * sizeof(T));
    return *this;
  }

  // Element access
  T& operator[](unsigned long i) { return data[i]; }
  const T& operator[](unsigned long i) const { return data[i]; }

  T& operator()(unsigned long i) {
    assert(i < n);
    return data[i];
  }

  const T& operator()(unsigned long i) const {
    assert(i < n);
    return data[i];
  }

  iterator begin() { return data; }
  const_iterator begin() const { return data; }

  iterator end() { return data + n; }
  const_iterator end() const { return data + n; }

  bool empty() const { return n == 0; }
  unsigned long capacity() const { return max_n; }

  void resize(unsigned long new_size) {
    if (new_size > max_n) reserve(new_size);
    n = new_size;
  }

  void reserve(unsigned long r) {
#ifndef NO_THROWS
    if (r > ULONG_MAX / sizeof(T)) throw std::bad_alloc();
#endif
    T* resized = static_cast<T*>(std::realloc(data, r * sizeof(T)));
#ifndef NO_THROWS
    if (!resized) throw std::bad_alloc();
#endif
    data = resized;
    max_n = r;
  }

  void set_zero() {
    if (data) std::memset(data, 0, n * sizeof(T));
  }

  void clear() {
    std::free(data);
    data = nullptr;
    n = max_n = 0;
  }
};

// WrapArray1<T>: Lightweight non-owning 1D array view
template <typename T>
struct WrapArray1 {
  using iterator = T*;
  using const_iterator = const T*;

  unsigned long n = 0;
  unsigned long max_n = 0;
  T* data = nullptr;

  WrapArray1() = default;

  WrapArray1(unsigned long n_, T* data_) : n(n_), max_n(n_), data(data_) { assert(data || max_n == 0); }

  WrapArray1(Array1<T>& a) : n(a.n), max_n(a.max_n), data(a.data) {}

  WrapArray1(std::vector<T>& a) : n(a.size()), max_n(a.capacity()), data(a.data()) {}

  T& operator[](unsigned long i) { return data[i]; }
  const T& operator[](unsigned long i) const { return data[i]; }

  T& operator()(unsigned long i) {
    assert(i < n);
    return data[i];
  }

  const T& operator()(unsigned long i) const {
    assert(i < n);
    return data[i];
  }

  iterator begin() { return data; }
  const_iterator begin() const { return data; }

  iterator end() { return data + n; }
  const_iterator end() const { return data + n; }

  void clear() { n = 0; }
  bool empty() const { return n == 0; }

  unsigned long size() const { return n; }
};