#pragma once

#include <vector>
#include "math_defs.h"

struct Particle;
class FluidSim;

/// Spatial hashing grid used to accelerate neighbor lookups
class sorter {
 public:
  /// Constructor sets up grid dimensions
  sorter(int ni_, int nj_);

  /// Destructor
  ~sorter();

  /// Assign particles to their corresponding cells based on position
  void sort(const std::vector<Particle>& particles, const Vector2s& origin, scalar dx);

  /// Visit all neighboring cells around (i, j) within given range offsets
  template <typename Callable>
  void getNeigboringParticles_cell(int i, int j, int wl, int wh, int hl, int hh, Callable func) {
    for (int sx = i + wl; sx <= i + wh; ++sx) {
      for (int sy = j + hl; sy <= j + hh; ++sy) {
        if (sx < 0 || sx >= ni || sy < 0 || sy >= nj) continue;
        int index = sy * ni + sx;
        func(cells[index]);
      }
    }
  }

  /// Return particle count in cell (i, j)
  int getNumParticleAt(int i, int j);

  std::vector<std::vector<const Particle*>> cells;  ///< 2D grid of particle pointers
  int ni;                                           ///< grid width
  int nj;                                           ///< grid height
};