#pragma once

#include <cmath>
#include "blas_wrapper.h"
#include "sparse_matrix.h"

namespace robertbridson {

//==============================================================
// A compressed sparse column format for a lower triangular matrix,
// separating out the diagonal entries for optimization purposes.
template <typename T>
struct SparseColumnLowerFactor {
  unsigned int n;
  std::vector<T> invdiag;
  std::vector<T> value;
  std::vector<unsigned int> rowindex;
  std::vector<unsigned int> colstart;
  std::vector<T> adiag;

  explicit SparseColumnLowerFactor(unsigned int n_ = 0) : n(n_), invdiag(n_), colstart(n_ + 1), adiag(n_) {}

  void resize(unsigned int new_size) {
    n = new_size;
    invdiag.resize(n);
    colstart.resize(n + 1);
    adiag.resize(n);
  }
};

//==============================================================
// Level-0 Incomplete Cholesky factorization with modification.
// The parameter `modification_parameter` controls the degree of modification,
// and `min_diagonal_ratio` is used to stabilize the diagonal.
template <typename T>
void factor_modified_incomplete_cholesky0(const SparseMatrix<T> &matrix, SparseColumnLowerFactor<T> &factor, T modification_parameter = 0.97,
                                          T min_diagonal_ratio = 0.25) {
  factor.resize(matrix.n);
  zero(factor.invdiag);
  zero(factor.adiag);
  factor.value.clear();
  factor.rowindex.clear();

  for (unsigned int i = 0; i < matrix.n; ++i) {
    factor.colstart[i] = static_cast<unsigned int>(factor.rowindex.size());
    for (unsigned int j = 0; j < matrix.index[i].size(); ++j) {
      unsigned int idx = matrix.index[i][j];
      T val = matrix.value[i][j];
      if (idx > i) {
        factor.rowindex.push_back(idx);
        factor.value.push_back(val);
      } else if (idx == i) {
        factor.invdiag[i] = val;
        factor.adiag[i] = val;
      }
    }
  }
  factor.colstart[matrix.n] = static_cast<unsigned int>(factor.rowindex.size());

  for (unsigned int k = 0; k < matrix.n; ++k) {
    if (factor.adiag[k] == 0) continue;

    if (factor.invdiag[k] < min_diagonal_ratio * factor.adiag[k])
      factor.invdiag[k] = 1 / std::sqrt(factor.adiag[k]);
    else
      factor.invdiag[k] = 1 / std::sqrt(factor.invdiag[k]);

    for (unsigned int p = factor.colstart[k]; p < factor.colstart[k + 1]; ++p) factor.value[p] *= factor.invdiag[k];

    for (unsigned int p = factor.colstart[k]; p < factor.colstart[k + 1]; ++p) {
      unsigned int j = factor.rowindex[p];
      T multiplier = factor.value[p];
      T missing = 0;
      unsigned int a = factor.colstart[k];
      unsigned int b = 0;

      while (a < factor.colstart[k + 1] && factor.rowindex[a] < j) {
        while (b < matrix.index[j].size()) {
          if (matrix.index[j][b] < factor.rowindex[a])
            ++b;
          else if (matrix.index[j][b] == factor.rowindex[a])
            break;
          else {
            missing += factor.value[a];
            break;
          }
        }
        ++a;
      }

      if (a < factor.colstart[k + 1] && factor.rowindex[a] == j) factor.invdiag[j] -= multiplier * factor.value[a];

      ++a;
      b = factor.colstart[j];

      while (a < factor.colstart[k + 1] && b < factor.colstart[j + 1]) {
        if (factor.rowindex[b] < factor.rowindex[a])
          ++b;
        else if (factor.rowindex[b] == factor.rowindex[a]) {
          factor.value[b] -= multiplier * factor.value[a];
          ++a;
          ++b;
        } else {
          missing += factor.value[a];
          ++a;
        }
      }

      while (a < factor.colstart[k + 1]) {
        missing += factor.value[a];
        ++a;
      }

      factor.invdiag[j] -= modification_parameter * multiplier * missing;
    }
  }
}

//==============================================================
// Solving L * x = rhs for lower triangular matrix L.
template <typename T>
void solve_lower(const SparseColumnLowerFactor<T> &factor, const std::vector<T> &rhs, std::vector<T> &result) {
  assert(factor.n == rhs.size() && factor.n == result.size());
  result = rhs;
  for (unsigned int i = 0; i < factor.n; ++i) {
    result[i] *= factor.invdiag[i];
    for (unsigned int j = factor.colstart[i]; j < factor.colstart[i + 1]; ++j) result[factor.rowindex[j]] -= factor.value[j] * result[i];
  }
}

// Solving L^T * x = rhs for upper triangular matrix from L^T.
template <typename T>
void solve_lower_transpose_in_place(const SparseColumnLowerFactor<T> &factor, std::vector<T> &x) {
  assert(factor.n == x.size() && factor.n > 0);
  for (int i = static_cast<int>(factor.n) - 1; i >= 0; --i) {
    for (unsigned int j = factor.colstart[i]; j < factor.colstart[i + 1]; ++j) x[i] -= factor.value[j] * x[factor.rowindex[j]];
    x[i] *= factor.invdiag[i];
  }
}

//==============================================================
// Preconditioned Conjugate Gradient solver with IC(0) preconditioner.
template <typename T>
struct PCGSolver {
  PCGSolver() { set_solver_parameters(1e-4, 100, 0.97, 0.25); }

  void set_solver_parameters(T tol, int max_iter, T mic_param = 0.97, T min_diag_ratio = 0.25) {
    tolerance_factor = std::max(tol, static_cast<T>(1e-30));
    max_iterations = max_iter;
    modified_incomplete_cholesky_parameter = mic_param;
    min_diagonal_ratio = min_diag_ratio;
  }

  bool solve(const SparseMatrix<T> &matrix, const std::vector<T> &rhs, std::vector<T> &result, T &residual_out, int &iterations_out) {
    const unsigned int n = matrix.n;
    if (m.size() != n) {
      m.resize(n);
      s.resize(n);
      z.resize(n);
      r.resize(n);
    }

    zero(result);
    r = rhs;
    residual_out = BLAS::abs_max(r);
    if (residual_out == 0) {
      iterations_out = 0;
      return true;
    }

    T tol = tolerance_factor * residual_out;
    form_preconditioner(matrix);
    apply_preconditioner(r, z);
    T rho = BLAS::dot(z, r);
    if (rho == 0 || std::isnan(rho)) {
      iterations_out = 0;
      return false;
    }

    s = z;
    fixed_matrix.construct_from_matrix(matrix);

    for (int iter = 0; iter < max_iterations; ++iter) {
      multiply(fixed_matrix, s, z);
      T alpha = rho / BLAS::dot(s, z);
      BLAS::add_scaled(alpha, s, result);
      BLAS::add_scaled(-alpha, z, r);

      residual_out = BLAS::abs_max(r);
      if (residual_out <= tol) {
        iterations_out = iter + 1;
        return true;
      }

      apply_preconditioner(r, z);
      T rho_new = BLAS::dot(z, r);
      T beta = rho_new / rho;
      BLAS::add_scaled(beta, s, z);
      s.swap(z);
      rho = rho_new;
    }

    iterations_out = max_iterations;
    return false;
  }

 protected:
  SparseColumnLowerFactor<T> ic_factor;
  std::vector<T> m, z, s, r;
  FixedSparseMatrix<T> fixed_matrix;

  T tolerance_factor;
  int max_iterations;
  T modified_incomplete_cholesky_parameter;
  T min_diagonal_ratio;

  void form_preconditioner(const SparseMatrix<T> &matrix) {
    factor_modified_incomplete_cholesky0(matrix, ic_factor, modified_incomplete_cholesky_parameter, min_diagonal_ratio);
  }

  void apply_preconditioner(const std::vector<T> &x, std::vector<T> &result) {
    solve_lower(ic_factor, x, result);
    solve_lower_transpose_in_place(ic_factor, result);
  }
};

}  