#pragma once

#include "math_defs.h"

namespace kernel {

/// Smooth cubic kernel (1 - r^2/h^2)^3 for r^2 < h^2
inline scalar smooth_kernel(const scalar& r2, const scalar& h) {
  const scalar q = r2 / (h * h);
  const scalar v = 1.0 - q;
  return v > 0.0 ? v * v * v : 0.0;
}

/// Laplacian of smooth kernel, linear decay for 0 <= r/h <= 1
inline scalar smooth_kernel_laplacian(const scalar& r2, const scalar& h) {
  const scalar normalized = sqrt(r2 / (h * h));
  return (normalized > 1.0) ? 0.0 : (1.0 - normalized);
}

/// Sharp kernel with inverse-square behavior capped at small r
inline scalar sharp_kernel(const scalar& r2, const scalar& h) {
  const scalar denom = std::max(r2, 1.0e-5f);
  const scalar term = (h * h) / denom - 1.0f;
  return term > 0.0f ? term : 0.0f;
}

/// Separable linear kernel in 2D: (1 - |dx/h|)(1 - |dy/h|)
inline scalar linear_kernel(const Vector2s& d, const scalar& h) {
  const scalar wx = 1.0 - std::abs(d(0) / h);
  const scalar wy = 1.0 - std::abs(d(1) / h);
  return (wx > 0.0 && wy > 0.0) ? wx * wy : 0.0;
}

/// 1D quadratic B-spline kernel
inline scalar quadratic_kernel_1d(const scalar& d, const scalar& h) {
  const scalar r = std::abs(d) / h;
  if (r < 0.5) {
    return 0.75 - r * r;
  } else if (r < 1.5) {
    const scalar t = 1.5 - r;
    return 0.5 * t * t;
  } else {
    return 0.0;
  }
}

/// 2D separable quadratic kernel: product of 1D kernels
inline scalar quadratic_kernel(const Vector2s& d, const scalar& h) { return quadratic_kernel_1d(d(0), h) * quadratic_kernel_1d(d(1), h); }

}  // namespace kernel