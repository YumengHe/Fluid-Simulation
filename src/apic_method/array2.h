#pragma once

#include <algorithm>
#include <cassert>
#include <vector>
#include "array1.h"


// 2D Array container based on a flat 1D array (row-major order).
// The default backend is std::vector<T>, but customizable (e.g., Array1<T>).
template <typename T, typename ArrayT = std::vector<T>>
struct Array2 {
  // Type aliases for compatibility
  using iterator = typename ArrayT::iterator;
  using const_iterator = typename ArrayT::const_iterator;
  using size_type = typename ArrayT::size_type;
  using reference = T&;
  using const_reference = const T&;
  using reverse_iterator = typename ArrayT::reverse_iterator;
  using const_reverse_iterator = typename ArrayT::const_reverse_iterator;

  int ni = 0;  // Number of columns
  int nj = 0;  // Number of rows
  ArrayT a;    // Flattened buffer (size = ni * nj)

  // Default constructor
  Array2() = default;

  // Element read access
  const T& operator()(int i, int j) const {
    assert(i >= 0 && j >= 0 && i < ni && j < nj);
    size_type idx = static_cast<size_type>(i + ni * j);
    return a[idx];
  }

  // Element write access
  T& operator()(int i, int j) {
    assert(i >= 0 && j >= 0 && i < ni && j < nj);
    size_type idx = static_cast<size_type>(i + ni * j);
    return a[idx];
  }

  // Iteration support
  inline iterator begin() { return a.begin(); }
  inline const_iterator begin() const { return a.begin(); }
  inline iterator end() { return a.end(); }
  inline const_iterator end() const { return a.end(); }

  // Utility functions
  inline bool empty() const { return a.empty(); }
  inline size_type size() const { return a.size(); }

  // Memory reservation
  void reserve(int reserve_ni, int reserve_nj) {
    const int total = reserve_ni * reserve_nj;
    a.reserve(total);
  }

  // Resize internal grid and update dimensions
  void resize(int new_ni, int new_nj) {
    assert(new_ni >= 0 && new_nj >= 0);
    ni = new_ni;
    nj = new_nj;
    a.resize(static_cast<size_type>(ni * nj));
  }

  // Set all entries to zero (assumes T is POD or supports memset-compatible zeroing)
  void set_zero() { a.set_zero(); }

  // Trim storage if backend supports it
  void trim() { a.trim(); }
};

// Frequently used 2D array types
using Array2s = Array2<scalar, Array1<scalar>>;
using Array2c = Array2<char, Array1<char>>;
using Array2i = Array2<int, Array1<int>>;

// Non-owning (wrapped) scalar grid
using WrapArray2s = Array2<scalar, WrapArray1<scalar>>;