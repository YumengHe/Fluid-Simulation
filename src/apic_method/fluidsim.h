// Copyright 2016 Raymond Yun Fei, Christopher Batty, Robert Bridson
//
// Licensed under the Apache License,
// Version 2.0(the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef APIC2D_FLUIDSIM_H_
#define APIC2D_FLUIDSIM_H_

#include <chrono>
#include <memory>
#include <vector>

#include "array2.h"
#include "math_defs.h"
#include "sparse_matrix.h"
#include "pcg_solver.h"

class sorter;

struct Particle {
  Particle(const Vector2s& x, const Vector2s& v, const scalar& radii, const scalar& density);
  Particle();
  Particle(const Particle&);

  Vector2s x_;
  Vector2s v_;
  Matrix2s c_;

  Vector2s buf0_;

  scalar radii_;
  scalar mass_;
  // scalar logJ_;
};

class FluidSim {
 public:
  using Clock = std::chrono::high_resolution_clock;
  using TimePoint = std::chrono::time_point<Clock>;
  using NeighborParticlesType = std::vector<const Particle*>;

  virtual ~FluidSim();

  scalar rho_;

  enum INTEGRATOR_TYPE {
    IT_APIC,

    IT_COUNT
  };

  enum VELOCITY_ORDER {
    VO_EULER,
    VO_RA2,
    VO_RK3,
    VO_RK4,

    VO_COUNT
  };

  enum INTERPOLATION_ORDER {
    IO_LINEAR,
    IO_QUADRATIC,

    IO_COUNT
  };

  enum BOUNDARY_TYPE {
    BT_BOX,

    BT_COUNT
  };

  struct Boundary {
    Boundary(const Vector2s& center, const Vector2s& parameter, BOUNDARY_TYPE type, bool inside);

    Boundary(Boundary* op0_, Boundary* op1_, BOUNDARY_TYPE type_);

    Vector2s center_;
    Vector2s parameter_;

    Boundary* op0_;
    Boundary* op1_;

    BOUNDARY_TYPE type_;
    scalar sign_;
  };


  void initialize(const Vector2s& origin, scalar width, int ni, int nj, scalar rho, bool draw_particles = true, bool print_timing = false);
  void advance(scalar dt);
  void update_boundary();
  void init_random_particles();
  void render();
  void render_boundaries(const Boundary& b);
  scalar compute_cfl() const;
  scalar solid_distance(const Vector2s& pos) const;
  scalar solid_distance(const Vector2s& pos, const Boundary& b) const;

  Vector2s get_velocity_and_affine_matrix_with_order(const Vector2s& position, scalar dt, FluidSim::VELOCITY_ORDER v_order,
                                                     FluidSim::INTERPOLATION_ORDER i_order, Matrix2s* affine_matrix);



  /*! Linear interpolation kernels */
  Vector2s get_velocity(const Vector2s& position);
  Matrix2s get_affine_matrix(const Vector2s& position);
  Vector2s get_saved_velocity(const Vector2s& position);
  Matrix2s get_saved_affine_matrix(const Vector2s& position);

  /*! Add particle to the system */
  void add_particle(const Particle& position);

  /*! P2G scheme */
  void map_p2g();
  void map_p2g_linear();

  /*! FLIP schemes */
  void map_g2p_flip_general(float dt, const scalar affine_stretching_ratio,
                            const scalar affine_rotational_ratio);

  void save_velocity();

  void relaxation(scalar dt);
  void set_root_boundary(const Boundary& b) { root_boundary_ = std::make_unique<Boundary>(b); }
  const std::vector<Particle>& get_particles() const { return particles_; }

  const Vector2s& get_origin() const { return origin_; }

 protected:
  inline scalar box_distance(const Vector2s& position, const Vector2s& centre, const Vector2s& expand) const {
    scalar dx = fabs(position[0] - centre[0]) - expand[0];
    scalar dy = fabs(position[1] - centre[1]) - expand[1];
    scalar dax = max(dx, 0.0f);
    scalar day = max(dy, 0.0f);
    return min(max(dx, dy), 0.0f) + sqrt(dax * dax + day * day);
  }

  inline void tick() {
    if (print_timing_) last_time_point_ = Clock::now();
  }

  inline void tock(const std::string& info) {
    if (print_timing_)
      std::cout << "[" << info << "] "
                << static_cast<scalar>(std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - last_time_point_).count()) * 1e-6 << " ms"
                << std::endl;
  }

  // tracer particle operations


  // fluid velocity operations
  void add_force(scalar dt);

  void compute_weights();
  void solve_pressure(scalar dt);
  void compute_liquid_distance();

  void constrain_velocity();

 private:
  /*! Boundaries */
  std::unique_ptr<Boundary> root_boundary_;

  /*! Grid Origin */
  Vector2s origin_;

  /*! Grid dimensions */
  int ni_, nj_;
  scalar dx_;

  /*! Fluid velocity */
  Array2s u_, v_;
  Array2s temp_u_, temp_v_;
  Array2s saved_u_, saved_v_;

  /*! Tracer particles */
  std::vector<Particle> particles_;

  /*! Static geometry representation */
  Array2s nodal_solid_phi_;
  Array2s liquid_phi_;
  Array2s u_weights_, v_weights_;

  /*! Data arrays for extrapolation */
  Array2c valid_, old_valid_;
  Array2c u_valid_, v_valid_;

  sorter* m_sorter_;

  /*! Solver data */
  robertbridson::PCGSolver<scalar> solver_;
  robertbridson::SparseMatrix<scalar> matrix_;
  std::vector<scalar> rhs_;
  std::vector<scalar> pressure_;

  bool draw_particles_;
  bool print_timing_;

  TimePoint last_time_point_;

};

#endif  // APIC2D_FLUIDSIM_H_
