#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>
#include <Eigen/StdVector>

// Scalar type used for floating-point operations
using scalar = float;

// Struct representing an (index, value) pair, with comparator
struct int_scalar {
  int i;
  scalar v;

  inline bool operator()(const int_scalar& lhs, const int_scalar& rhs) { return lhs.v < rhs.v; }
};

// === Commonly used vector and matrix types ===

// 2D float vector (column) and its integer counterpart
using Vector2s = Eigen::Matrix<scalar, 2, 1>;
using Vector2i = Eigen::Matrix<int, 2, 1>;

// Transposed 2D float vector (row)
using Vector2sT = Eigen::Matrix<scalar, 1, 2>;

// 2x2 float matrix
using Matrix2s = Eigen::Matrix<scalar, 2, 2>;

// Dynamically-sized float vector
using VectorXs = Eigen::Matrix<scalar, Eigen::Dynamic, 1>;

// === Sparse matrix utilities ===

// Sparse float matrix
using SparseXs = Eigen::SparseMatrix<scalar>;

// Triplet used to construct sparse matrices
using Triplets = Eigen::Triplet<scalar>;
using TripletXs = std::vector<Triplets>;