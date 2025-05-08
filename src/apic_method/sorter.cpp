#include "sorter.h"
#include "fluidsim.h"

using std::max;
using std::min;
using std::vector;

// Constructor initializes cell grid
sorter::sorter(int ni_, int nj_) : ni(ni_), nj(nj_) { cells.resize(ni * nj); }

sorter::~sorter() = default;

/// Assign each particle to its grid cell for spatial lookup
void sorter::sort(const vector<Particle>& particles, const Vector2s& origin, scalar dx) {
  // Clear previous cell assignments
  for (int y = 0; y < nj; ++y) {
    for (int x = 0; x < ni; ++x) {
      int idx = y * ni + x;
      cells[idx].clear();
    }
  }

  // Insert each particle into its corresponding cell
  const int particle_count = static_cast<int>(particles.size());
  for (int n = 0; n < particle_count; ++n) {
    const Particle& p = particles[n];

    scalar px = p.x_(0);
    scalar py = p.x_(1);

    int cell_x = static_cast<int>((px - origin(0)) / dx);
    int cell_y = static_cast<int>((py - origin(1)) / dx);

    int clamped_x = max(0, min(ni - 1, cell_x));
    int clamped_y = max(0, min(nj - 1, cell_y));

    int flat_index = clamped_y * ni + clamped_x;
    cells[flat_index].push_back(&p);
  }
}
