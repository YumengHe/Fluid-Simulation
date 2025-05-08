#pragma once

#include <Eigen/Core>
#include <cmath>
#include <vector>

namespace robertbridson {
namespace BLAS {

/// Compute dot product of two float vectors
inline float dot(const std::vector<float>& x, const std::vector<float>& y) {
  const auto vecX = Eigen::Map<const Eigen::VectorXf>(x.data(), static_cast<int>(x.size()));
  const auto vecY = Eigen::Map<const Eigen::VectorXf>(y.data(), static_cast<int>(y.size()));
  return vecX.dot(vecY);
}

/// Find index of element with maximum absolute value (infinity-norm location)
inline int index_abs_max(const std::vector<float>& x) {
  const auto vec = Eigen::Map<const Eigen::VectorXf>(x.data(), static_cast<int>(x.size()));
  int index = 0;
  vec.maxCoeff(&index);  // Eigen will fill index
  return index;
}

/// Compute the maximum absolute value (infinity-norm value)
inline float abs_max(const std::vector<float>& x) {
  const int idx = index_abs_max(x);
  return std::fabs(x[idx]);
}

/// Perform scaled vector addition: y = alpha * x + y
inline void add_scaled(float alpha, const std::vector<float>& x, std::vector<float>& y) {
  Eigen::Map<Eigen::VectorXf> vecY(y.data(), static_cast<int>(y.size()));
  const auto vecX = Eigen::Map<const Eigen::VectorXf>(x.data(), static_cast<int>(x.size()));
  vecY += alpha * vecX;
}

}  // namespace BLAS
}  // namespace robertbridson