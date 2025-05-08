#pragma once

#include <cassert>
#include <vector>

#include "util.h"

namespace robertbridson {

/// A dynamic Compressed Sparse Row (CSR) matrix implementation
template <typename T>
class SparseMatrix {
 public:
  unsigned int n;                                // number of rows/columns
  std::vector<std::vector<unsigned int>> index;  // row-wise column indices
  std::vector<std::vector<T>> value;             // row-wise values

  explicit SparseMatrix(unsigned int dim = 0, unsigned int estimate_per_row = 7) : n(dim), index(dim), value(dim) {
    for (unsigned int i = 0; i < dim; ++i) {
      index[i].reserve(estimate_per_row);
      value[i].reserve(estimate_per_row);
    }
  }

  void clear() {
    n = 0;
    index.clear();
    value.clear();
  }

  void zero() {
    for (unsigned int i = 0; i < n; ++i) {
      index[i].clear();
      value[i].clear();
    }
  }

  void resize(int new_n) {
    n = new_n;
    index.resize(n);
    value.resize(n);
  }

  T operator()(unsigned int i, unsigned int j) const {
    const auto& row = index[i];
    const auto& val = value[i];
    for (size_t k = 0; k < row.size(); ++k) {
      if (row[k] == j) return val[k];
      if (row[k] > j) break;
    }
    return T(0);
  }

  void set_element(unsigned int i, unsigned int j, T val) {
    auto& row_idx = index[i];
    auto& row_val = value[i];

    for (size_t k = 0; k < row_idx.size(); ++k) {
      if (row_idx[k] == j) {
        row_val[k] = val;
        return;
      } else if (row_idx[k] > j) {
        insert(row_idx, k, j);
        insert(row_val, k, val);
        return;
      }
    }

    row_idx.push_back(j);
    row_val.push_back(val);
  }

  void add_to_element(unsigned int i, unsigned int j, T increment) {
    auto& row_idx = index[i];
    auto& row_val = value[i];

    for (size_t k = 0; k < row_idx.size(); ++k) {
      if (row_idx[k] == j) {
        row_val[k] += increment;
        return;
      } else if (row_idx[k] > j) {
        insert(row_idx, k, j);
        insert(row_val, k, increment);
        return;
      }
    }

    row_idx.push_back(j);
    row_val.push_back(increment);
  }
};

using SparseMatrixd = SparseMatrix<double>;

/// Perform matrix-vector product: result = A * x
template <typename T>
void multiply(const SparseMatrix<T>& A, const std::vector<T>& x, std::vector<T>& result) {
  assert(A.n == x.size());
  result.assign(A.n, T(0));
  for (unsigned int i = 0; i < A.n; ++i) {
    for (size_t j = 0; j < A.index[i].size(); ++j) {
      result[i] += A.value[i][j] * x[A.index[i][j]];
    }
  }
}

/// Perform result -= A * x
template <typename T>
void multiply_and_subtract(const SparseMatrix<T>& A, const std::vector<T>& x, std::vector<T>& result) {
  assert(A.n == x.size());
  for (unsigned int i = 0; i < A.n; ++i) {
    for (size_t j = 0; j < A.index[i].size(); ++j) {
      result[i] -= A.value[i][j] * x[A.index[i][j]];
    }
  }
}

/// Fixed compressed sparse row format (CSR), more efficient for mat-vec
template <typename T>
class FixedSparseMatrix {
 public:
  unsigned int n;                      // number of rows
  std::vector<T> value;                // non-zero values
  std::vector<unsigned int> colindex;  // corresponding column indices
  std::vector<unsigned int> rowstart;  // row start positions

  explicit FixedSparseMatrix(unsigned int dim = 0) : n(dim), value(), colindex(), rowstart(dim + 1) {}

  void resize(int new_n) {
    n = new_n;
    rowstart.resize(n + 1);
  }

  void construct_from_matrix(const SparseMatrix<T>& source) {
    resize(source.n);
    rowstart[0] = 0;

    for (unsigned int i = 0; i < n; ++i) {
      rowstart[i + 1] = rowstart[i] + source.index[i].size();
    }

    value.resize(rowstart[n]);
    colindex.resize(rowstart[n]);

    unsigned int cursor = 0;
    for (unsigned int i = 0; i < n; ++i) {
      for (size_t k = 0; k < source.index[i].size(); ++k) {
        value[cursor] = source.value[i][k];
        colindex[cursor] = source.index[i][k];
        ++cursor;
      }
    }
  }
};

using FixedSparseMatrixd = FixedSparseMatrix<double>;

/// Multiply fixed sparse matrix with a vector: result = A * x
template <typename T>
void multiply(const FixedSparseMatrix<T>& A, const std::vector<T>& x, std::vector<T>& result) {
  assert(A.n == x.size());
  result.assign(A.n, T(0));

  for (unsigned int i = 0; i < A.n; ++i) {
    for (unsigned int j = A.rowstart[i]; j < A.rowstart[i + 1]; ++j) {
      result[i] += A.value[j] * x[A.colindex[j]];
    }
  }
}

}  // namespace robertbridson