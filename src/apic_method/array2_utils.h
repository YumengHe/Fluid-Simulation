#pragma once

#include "array2.h"
#include "util.h"

// Interpolate a value from a 2D scalar grid using bilinear interpolation
template <typename S, typename T>
T interpolate_value(const Eigen::Matrix<S, 2, 1>& pos, const Array2<T, Array1<T>>& grid) {
  int ix, iy;
  S wx, wy;

  get_barycentric(pos[0], ix, wx, 0, grid.ni);
  get_barycentric(pos[1], iy, wy, 0, grid.nj);

  const T& v00 = grid(ix, iy);
  const T& v10 = grid(ix + 1, iy);
  const T& v01 = grid(ix, iy + 1);
  const T& v11 = grid(ix + 1, iy + 1);

  const auto lerp_x0 = (1 - wx) * v00 + wx * v10;
  const auto lerp_x1 = (1 - wx) * v01 + wx * v11;

  return (1 - wy) * lerp_x0 + wy * lerp_x1;
}

// Compute affine-interpolated gradient from a 2D grid
template <typename T>
Eigen::Matrix<T, 2, 1> affine_interpolate_value(const Eigen::Matrix<T, 2, 1>& pos, const Array2<T, Array1<T>>& grid) {
  int ix, iy;
  T wx, wy;

  get_barycentric(pos[0], ix, wx, 0, grid.ni);
  get_barycentric(pos[1], iy, wy, 0, grid.nj);

  const T& v00 = grid(ix, iy);
  const T& v10 = grid(ix + 1, iy);
  const T& v01 = grid(ix, iy + 1);
  const T& v11 = grid(ix + 1, iy + 1);

  T dx0 = v10 - v00;
  T dx1 = v11 - v01;
  T dy0 = v01 - v00;
  T dy1 = v11 - v10;

  T dx = dx0 * (1 - wy) + dx1 * wy;
  T dy = dy0 * (1 - wx) + dy1 * wx;

  return Eigen::Matrix<T, 2, 1>(dx, dy);
}

// Compute interpolated value and gradient from a scalar grid
template <typename S, typename T>
T interpolate_gradient(Eigen::Matrix<T, 2, 1>& grad, const Eigen::Matrix<S, 2, 1>& pos, const Array2<T, Array1<T>>& grid) {
  int ix, iy;
  S wx, wy;

  get_barycentric(pos[0], ix, wx, 0, grid.ni);
  get_barycentric(pos[1], iy, wy, 0, grid.nj);

  const T& v00 = grid(ix, iy);
  const T& v10 = grid(ix + 1, iy);
  const T& v01 = grid(ix, iy + 1);
  const T& v11 = grid(ix + 1, iy + 1);

  T gx0 = v10 - v00;
  T gx1 = v11 - v01;
  T gy0 = v01 - v00;
  T gy1 = v11 - v10;

  grad[0] = gx0 * (1 - wy) + gx1 * wy;
  grad[1] = gy0 * (1 - wx) + gy1 * wx;

  const auto interp_x0 = (1 - wx) * v00 + wx * v10;
  const auto interp_x1 = (1 - wx) * v01 + wx * v11;

  return (1 - wy) * interp_x0 + wy * interp_x1;
}