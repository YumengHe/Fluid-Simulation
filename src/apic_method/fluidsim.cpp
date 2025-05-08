// ==== Headers and Configuration ====
#include "fluidsim.h"

#include "array2_utils.h"
#include "kernel.h"
#include "sorter.h"
#include "sparse_matrix.h"


#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

//  Simulation Scheme Settings 
const FluidSim::INTEGRATOR_TYPE integration_scheme = FluidSim::IT_APIC;
const FluidSim::VELOCITY_ORDER velocity_order = FluidSim::VO_EULER;
const FluidSim::INTERPOLATION_ORDER interpolation_order = FluidSim::IO_LINEAR;

//  Particle Behavior Settings 
const scalar lagrangian_ratio = 0.97f;
const int particle_correction_step = 1;

//  Function Prototypes 
scalar fraction_inside(scalar phi_left, scalar phi_right);
void extrapolate(Array2s& grid, Array2s& old_grid, const Array2s& grid_weight, 
                 const Array2s& grid_liquid_weight, Array2c& valid_, Array2c old_valid_,
                 const Vector2i& offset, int num_layers);

//  Destructor Implementation 
FluidSim::~FluidSim() {
  if (m_sorter_) {
    delete m_sorter_;
  }
}

void FluidSim::initialize(const Vector2s& origin, scalar width, int ni, int nj, scalar rho, bool draw_particles, bool print_timing) {
  // Store basic domain configuration
  origin_ = origin;
  ni_ = ni;
  nj_ = nj;
  dx_ = width / static_cast<scalar>(ni);
  rho_ = rho;
  draw_particles_ = draw_particles;
  print_timing_ = print_timing;

  // Resize velocity-related grid fields
  u_.resize(ni_ + 1, nj_);
  temp_u_.resize(ni_ + 1, nj_);
  u_weights_.resize(ni_ + 1, nj_);
  u_valid_.resize(ni_ + 1, nj_);
  saved_u_.resize(ni_ + 1, nj_);

  v_.resize(ni_, nj_ + 1);
  temp_v_.resize(ni_, nj_ + 1);
  v_weights_.resize(ni_, nj_ + 1);
  v_valid_.resize(ni_, nj_ + 1);
  saved_v_.resize(ni_, nj_ + 1);

  // Zero velocity field
  u_.set_zero();
  v_.set_zero();

  // Resize scalar and marker fields
  nodal_solid_phi_.resize(ni_ + 1, nj_ + 1);
  valid_.resize(ni_ + 1, nj_ + 1);
  old_valid_.resize(ni_ + 1, nj_ + 1);
  liquid_phi_.resize(ni_, nj_);

  // Allocate sorter for particle-neighbor queries
  m_sorter_ = new sorter(ni_, nj_);
}

// Update the signed distance field for all nodes in the scalar grid
void FluidSim::update_boundary() {
  parallel_for(0, nj_ + 1, [this](int j) {
    for (int i = 0; i <= ni_; ++i) {
      Vector2s pos_grid = Vector2s(i * dx_, j * dx_) + origin_;
      nodal_solid_phi_(i, j) = solid_distance(pos_grid);
    }
  });
}

// Perform particle relaxation using spring-like interactions with neighbors
void FluidSim::relaxation(scalar dt) {
  int particle_count = static_cast<int>(particles_.size());
  const scalar influence_radius = dx_ / std::sqrt(2.0) * 1.1f;
  int phase_offset = rand() % particle_correction_step;

  // Predict smoothed positions based on neighboring particle forces
  parallel_for(0, particle_count, [&](int index) {
    if (index % particle_correction_step != phase_offset) return;

    Particle& particle = particles_[index];
    Vector2s force_accum = Vector2s::Zero();

    int i = clamp(static_cast<int>((particle.x_(0) - origin_(0)) / dx_), 0, ni_);
    int j = clamp(static_cast<int>((particle.x_(1) - origin_(1)) / dx_), 0, nj_);

    m_sorter_->getNeigboringParticles_cell(i, j, -1, 1, -1, 1, [&](const NeighborParticlesType& neighbors) {
      for (const Particle* neighbor : neighbors) {
        if (&particle == neighbor) continue;
        scalar distance = (particle.x_ - neighbor->x_).norm();
        scalar kernel_weight = 50.0f * kernel::smooth_kernel(distance * distance, influence_radius);

        if (distance > 0.01f * influence_radius) {
          force_accum += kernel_weight * (particle.x_ - neighbor->x_) / distance * influence_radius;
        } else {
          force_accum(0) += 0.01f * influence_radius / dt * ((rand() & 0xFF) / 255.0f);
          force_accum(1) += 0.01f * influence_radius / dt * ((rand() & 0xFF) / 255.0f);
        }
      }
    });

    particle.buf0_ = particle.x_ + dt * force_accum;

    Vector2s grid_pos = (particle.buf0_ - origin_) / dx_;
    scalar phi_val = interpolate_value(grid_pos, nodal_solid_phi_);

    if (phi_val < 0) {
      Vector2s grad;
      interpolate_gradient(grad, grid_pos, nodal_solid_phi_);
      grad.normalize();
      particle.buf0_ -= phi_val * grad;
    }
  });

  // Commit predicted positions to actual particle state
  parallel_for(0, particle_count, [&](int index) {
    if (index % particle_correction_step != phase_offset) return;
    particles_[index].x_ = particles_[index].buf0_;
  });
}

/*!
  \brief  Executes one full simulation step: transfers data, applies forces, solves pressure, updates particles.
*/
void FluidSim::advance(scalar dt) {
  // Step 1: Sort particles spatially for neighbor search
  tick();
  m_sorter_->sort(particles_, origin_, dx_);
  tock("sort");

  // Step 2: Transfer particle velocities to the grid (P2G)
  tick();
  map_p2g();
  tock("p2g");

  // Step 3: Backup current grid velocities
  tick();
  save_velocity();
  tock("save velocity");

  // Step 4: Apply external forces like gravity
  tick();
  add_force(dt);
  tock("add force");

  // Step 5: Compute distance field to the fluid surface
  tick();
  compute_liquid_distance();
  tock("compute phi");

  // Step 6: Calculate face weights for fluid simulation
  tick();
  compute_weights();
  tock("compute weights");

  // Step 7: Solve for pressure to ensure incompressibility
  tick();
  solve_pressure(dt);
  tock("solve pressure");

  // Step 8: Extrapolate velocity field into invalid regions (faces with zero weight)
  tick();
  extrapolate(u_, temp_u_, u_weights_, liquid_phi_, valid_, old_valid_, Vector2i(-1, 0), 2);
  extrapolate(v_, temp_v_, v_weights_, liquid_phi_, valid_, old_valid_, Vector2i(0, -1), 2);
  tock("extrapolate");

  // Step 9: Enforce boundary conditions (e.g. solid wall velocity normal zero)
  tick();
  constrain_velocity();
  tock("constrain velocity");

  // Step 10: Smooth particle distribution using spring relaxation
  tick();
  relaxation(dt);
  tock("relaxation");

  // Step 11: Transfer velocities back from grid to particles (G2P), with affine APIC enhancement
  tick();
  map_g2p_flip_general(dt, 1.0, 1.0);
  tock("g2p");
}

// Store the current velocity field if Lagrangian blending is enabled
void FluidSim::save_velocity() {
  if (lagrangian_ratio > 0.0) {
    saved_u_ = u_;  // Backup horizontal component
    saved_v_ = v_;  // Backup vertical component
  }
}

// Apply external acceleration (e.g., gravity) to the vertical velocity field
void FluidSim::add_force(scalar dt) {
  for (int j = 0; j <= nj_; ++j) {
    for (int i = 0; i < ni_; ++i) {
      // Apply downward gravitational force to each vertical face
      v_(i, j) += -981.0f * dt;
    }
  }
}

// Estimate maximum stable timestep using CFL condition
scalar FluidSim::compute_cfl() const {
  scalar max_velocity_sq = 0.0;
  for (const Particle& p : particles_) {
    max_velocity_sq = std::max(max_velocity_sq, p.v_.squaredNorm());
  }
  return dx_ / std::sqrt(max_velocity_sq);
}

// Zero out the velocity component normal to solid boundaries in extrapolated regions
void FluidSim::constrain_velocity() {
  // Loop over the larger of the u and v grid dimensions
  parallel_for(0, std::max(u_.nj, v_.nj), [this](int j) {
    // Adjust u velocities
    if (j < u_.nj) {
      for (int i = 0; i < u_.ni; ++i) {
        if (u_weights_(i, j) == 0) {
          Vector2s grid_pos(i * dx_, (j + 0.5f) * dx_);
          Vector2s world_pos = grid_pos + origin_;
          Vector2s vel = get_velocity(world_pos);
          Vector2s normal_vec;
          interpolate_gradient(normal_vec, Vector2s(i, j + 0.5f), nodal_solid_phi_);
          if (normal_vec.norm() > 0.0f) normal_vec.normalize();
          scalar projected = vel.dot(normal_vec);
          u_(i, j) = vel[0] - projected * normal_vec[0];
        }
      }
    }

    // Adjust v velocities
    if (j < v_.nj) {
      for (int i = 0; i < v_.ni; ++i) {
        if (v_weights_(i, j) == 0) {
          Vector2s grid_pos((i + 0.5f) * dx_, j * dx_);
          Vector2s world_pos = grid_pos + origin_;
          Vector2s vel = get_velocity(world_pos);
          Vector2s normal_vec;
          interpolate_gradient(normal_vec, Vector2s(i + 0.5f, j), nodal_solid_phi_);
          if (normal_vec.norm() > 0.0f) normal_vec.normalize();
          scalar projected = vel.dot(normal_vec);
          v_(i, j) = vel[1] - projected * normal_vec[1];
        }
      }
    }
  });
}

// Build signed distance field for the fluid region based on nearby particles
void FluidSim::compute_liquid_distance() {
  const scalar min_radius = dx_ / std::sqrt(2.0f);

  parallel_for(0, static_cast<int>(nj_), [&](int j) {
    for (int i = 0; i < ni_; ++i) {
      Vector2s cell_center = Vector2s((i + 0.5f) * dx_, (j + 0.5f) * dx_) + origin_;

      scalar best_phi = dx_;  // initialize with max possible value

      // Check neighboring particles
      m_sorter_->getNeigboringParticles_cell(i, j, -1, 1, -1, 1, [&](const NeighborParticlesType& neighbors) {
        for (const Particle* p : neighbors) {
          scalar r = std::max(p->radii_, min_radius);
          scalar dist = (cell_center - p->x_).norm() - r;
          best_phi = std::min(best_phi, dist);
        }
      });

      // Apply solid boundary influence
      scalar solid_phi = solid_distance(cell_center);
      liquid_phi_(i, j) = std::min(best_phi, solid_phi);
    }
  });
}

// Compute velocity and optionally the affine C matrix at a given position, using specified evaluation orders
Vector2s FluidSim::get_velocity_and_affine_matrix_with_order(const Vector2s& position, scalar dt, FluidSim::VELOCITY_ORDER v_order,
                                                             FluidSim::INTERPOLATION_ORDER i_order, Matrix2s* affine_matrix) {
  auto velocity_func = &FluidSim::get_velocity;
  auto matrix_func = &FluidSim::get_affine_matrix;

  if (affine_matrix) {
    *affine_matrix = (this->*matrix_func)(position);
  }
  return (this->*velocity_func)(position);
}

// Sample velocity field from the MAC grid at a given point in space
Vector2s FluidSim::get_velocity(const Vector2s& position) {
  Vector2s local_pos = (position - origin_) / dx_;
  Vector2s u_pos = local_pos - Vector2s(0.0, 0.5);
  Vector2s v_pos = local_pos - Vector2s(0.5, 0.0);

  scalar u_interp = interpolate_value(u_pos, u_);
  scalar v_interp = interpolate_value(v_pos, v_);

  return Vector2s(u_interp, v_interp);
}

// Compute the affine velocity matrix (C matrix) for a position using grid interpolation
Matrix2s FluidSim::get_affine_matrix(const Vector2s& position) {
  Vector2s grid_pos = (position - origin_) / dx_;
  Vector2s u_loc = grid_pos - Vector2s(0.0, 0.5);
  Vector2s v_loc = grid_pos - Vector2s(0.5, 0.0);

  Matrix2s result;
  result.col(0) = affine_interpolate_value(u_loc, u_) / dx_;
  result.col(1) = affine_interpolate_value(v_loc, v_) / dx_;
  return result;
}

// Calculate how much of a segment is inside the fluid, based on signed distance values
scalar fraction_inside(scalar phi_a, scalar phi_b) {
  if (phi_a < 0 && phi_b < 0) return 1.0;
  if (phi_a < 0 && phi_b >= 0) return phi_a / (phi_a - phi_b);
  if (phi_a >= 0 && phi_b < 0) return phi_b / (phi_b - phi_a);
  return 0.0;
}

// Determine face weights based on solid boundaries for finite-volume treatment
void FluidSim::compute_weights() {
  parallel_for(0, std::max(u_weights_.nj, v_weights_.nj), [this](int j) {
    if (j < u_weights_.nj) {
      for (int i = 0; i < u_weights_.ni; ++i) {
        scalar phi1 = nodal_solid_phi_(i, j + 1);
        scalar phi2 = nodal_solid_phi_(i, j);
        u_weights_(i, j) = 1.0f - fraction_inside(phi1, phi2);
        u_weights_(i, j) = clamp(u_weights_(i, j), 0.0f, 1.0f);
      }
    }

    if (j < v_weights_.nj) {
      for (int i = 0; i < v_weights_.ni; ++i) {
        scalar phi1 = nodal_solid_phi_(i + 1, j);
        scalar phi2 = nodal_solid_phi_(i, j);
        v_weights_(i, j) = 1.0f - fraction_inside(phi1, phi2);
        v_weights_(i, j) = clamp(v_weights_(i, j), 0.0f, 1.0f);
      }
    }
  });
}

// Variational pressure projection solve for static geometry
void FluidSim::solve_pressure(scalar dt) {
  int N = ni_ * nj_;
  if (rhs_.size() != N) {
    rhs_.resize(N);
    pressure_.resize(N);
    matrix_.resize(N);
  }
  matrix_.zero();

  // Build linear system
  parallel_for(1, nj_ - 1, [&](int j) {
    for (int i = 1; i < ni_ - 1; ++i) {
      int idx = i + ni_ * j;
      rhs_[idx] = 0;
      pressure_[idx] = 0;
      float phi_c = liquid_phi_(i, j);
      if (phi_c < 0 && (u_weights_(i, j) > 0.0 || u_weights_(i + 1, j) > 0.0 || v_weights_(i, j) > 0.0 || v_weights_(i, j + 1) > 0.0)) {
        auto add_term = [&](float weight, float phi_n, int offset) {
          float coeff = weight * dt / sqr(dx_);
          if (phi_n < 0) {
            matrix_.add_to_element(idx, idx, coeff);
            matrix_.add_to_element(idx, idx + offset, -coeff);
          } else {
            float theta = std::max(0.01f, fraction_inside(phi_c, phi_n));
            matrix_.add_to_element(idx, idx, coeff / theta);
          }
        };

        add_term(u_weights_(i + 1, j), liquid_phi_(i + 1, j), +1);
        rhs_[idx] -= u_weights_(i + 1, j) * u_(i + 1, j) / dx_;

        add_term(u_weights_(i, j), liquid_phi_(i - 1, j), -1);
        rhs_[idx] += u_weights_(i, j) * u_(i, j) / dx_;

        add_term(v_weights_(i, j + 1), liquid_phi_(i, j + 1), +ni_);
        rhs_[idx] -= v_weights_(i, j + 1) * v_(i, j + 1) / dx_;

        add_term(v_weights_(i, j), liquid_phi_(i, j - 1), -ni_);
        rhs_[idx] += v_weights_(i, j) * v_(i, j) / dx_;
      }
    }
  });

  // Solve system
  scalar residual;
  int iters;
  bool ok = solver_.solve(matrix_, rhs_, pressure_, residual, iters);
  if (!ok) {
    std::cout << "WARNING: Pressure solve failed! residual = " << residual << ", iters = " << iters << std::endl;
  }

  // Apply pressure gradient to velocities
  parallel_for(0, std::max(u_.nj, v_.nj - 1), [&](int j) {
    if (j < u_.nj) {
      for (int i = 1; i < u_.ni - 1; ++i) {
        int idx = i + j * ni_;
        if (u_weights_(i, j) > 0) {
          if (liquid_phi_(i, j) < 0 || liquid_phi_(i - 1, j) < 0) {
            float theta = std::max(0.01f, fraction_inside(liquid_phi_(i - 1, j), liquid_phi_(i, j)));
            u_(i, j) -= dt * (pressure_[idx] - pressure_[idx - 1]) / dx_ / theta;
            u_valid_(i, j) = 1;
          } else {
            u_valid_(i, j) = 0;
          }
        } else {
          u_(i, j) = 0;
          u_valid_(i, j) = 0;
        }
      }
    }

    if (j >= 1 && j < v_.nj - 1) {
      for (int i = 0; i < v_.ni; ++i) {
        int idx = i + j * ni_;
        if (v_weights_(i, j) > 0) {
          if (liquid_phi_(i, j) < 0 || liquid_phi_(i, j - 1) < 0) {
            float theta = std::max(0.01f, fraction_inside(liquid_phi_(i, j - 1), liquid_phi_(i, j)));
            v_(i, j) -= dt * (pressure_[idx] - pressure_[idx - ni_]) / dx_ / theta;
            v_valid_(i, j) = 1;
          } else {
            v_valid_(i, j) = 0;
          }
        } else {
          v_(i, j) = 0;
          v_valid_(i, j) = 0;
        }
      }
    }
  });
}

scalar FluidSim::solid_distance(const Vector2s& pos, const Boundary& b) const {
  // Use box_distance multiplied by the sign to determine inside/outside
  return static_cast<scalar>(b.sign_) * box_distance(pos, b.center_, b.parameter_);
}

scalar FluidSim::solid_distance(const Vector2s& pos) const {
  // Apply default boundary if not specified
  return this->solid_distance(pos, *root_boundary_);
}

void FluidSim::init_random_particles() {
  int total_particles = ni_ * nj_;
  scalar radius = dx_ / std::sqrt(2.0);
  scalar spacing = dx_ * 20.0;

  for (int i = 0; i < ni_; ++i) {
    for (int j = 0; j < nj_; ++j) {
      for (int k = 0; k < 2; ++k) {
        scalar rand_x = static_cast<scalar>(rand()) / static_cast<scalar>(RAND_MAX);
        scalar rand_y = static_cast<scalar>(rand()) / static_cast<scalar>(RAND_MAX);
        scalar jitter_x = (rand_x * 2.0f - 1.0f + 0.5f) * dx_;
        scalar jitter_y = (rand_y * 2.0f - 1.0f + 0.5f) * dx_;
        Vector2s pos = origin_ + Vector2s(i * dx_, j * dx_) + Vector2s(jitter_x, jitter_y);

        if (solid_distance(pos) > spacing) {
          particles_.emplace_back(pos, Vector2s::Zero(), radius, rho_);
        }
      }
    }
  }
}

void FluidSim::map_p2g() {
  // Delegate to the linear implementation
  this->map_p2g_linear();
}

void FluidSim::map_p2g_linear() {
  // Project particle data onto the grid for both u and v velocity components
  parallel_for(0, nj_ + 1, [this](int j) {
    // Handle u-component: horizontal face-centered velocity
    if (j < nj_) {
      for (int i = 0; i <= ni_; ++i) {
        Vector2s grid_pos(i * dx_, (j + 0.5f) * dx_);
        Vector2s world_pos = grid_pos + origin_;

        scalar weighted_sum = 0.0;
        scalar total_weight = 0.0;

        m_sorter_->getNeigboringParticles_cell(i, j, -1, 0, -1, 1, [&](const NeighborParticlesType& neighbors) {
          for (const Particle* particle : neighbors) {
            Vector2s diff = world_pos - particle->x_;
            scalar weight = particle->mass_ * kernel::linear_kernel(diff, dx_);
            weighted_sum += weight * (particle->v_(0) + particle->c_.col(0).dot(diff));
            total_weight += weight;
          }
        });

        u_(i, j) = (total_weight > 0.0) ? weighted_sum / total_weight : 0.0;
      }
    }

    // Handle v-component: vertical face-centered velocity
    for (int i = 0; i < ni_; ++i) {
      Vector2s grid_pos((i + 0.5f) * dx_, j * dx_);
      Vector2s world_pos = grid_pos + origin_;

      scalar weighted_sum = 0.0;
      scalar total_weight = 0.0;

      m_sorter_->getNeigboringParticles_cell(i, j, -1, 1, -1, 0, [&](const NeighborParticlesType& neighbors) {
        for (const Particle* particle : neighbors) {
          Vector2s diff = world_pos - particle->x_;
          scalar weight = particle->mass_ * kernel::linear_kernel(diff, dx_);
          weighted_sum += weight * (particle->v_(1) + particle->c_.col(1).dot(diff));
          total_weight += weight;
        }
      });

      v_(i, j) = (total_weight > 0.0) ? weighted_sum / total_weight : 0.0;
    }
  });
}

void FluidSim::map_g2p_flip_general(float dt, const scalar affine_stretching_ratio, const scalar affine_rotational_ratio) {
  const bool use_affine = (affine_stretching_ratio > 0.0f || affine_rotational_ratio > 0.0f);

  // Update particles from grid (G2P) using affine FLIP variant
  parallel_for(0, static_cast<int>(particles_.size()), [&](int i) {
    Particle& particle = particles_[i];
    Matrix2s C_update = Matrix2s::Zero();

    // Retrieve updated velocity and affine matrix from grid
    Vector2s velocity_from_grid =
        get_velocity_and_affine_matrix_with_order(particle.x_, dt, velocity_order, interpolation_order, use_affine ? &C_update : nullptr);

    // Apply velocity and affine matrix back to the particle
    particle.v_ = velocity_from_grid;
    particle.c_ = C_update;
    particle.x_ += dt * velocity_from_grid;
  });
}

void FluidSim::render_boundaries(const Boundary& b) {
  Vector2s lower = b.center_ - b.parameter_;
  Vector2s upper = b.center_ + b.parameter_;

  // Draw rectangular boundary using line loop
  glColor3f(1.0f, 0.0f, 0.0f);  // red color
  glBegin(GL_LINE_LOOP);
  glVertex2f(lower(0), lower(1));
  glVertex2f(upper(0), lower(1));
  glVertex2f(upper(0), upper(1));
  glVertex2f(lower(0), upper(1));
  glEnd();
}

void FluidSim::render() {
  if (!draw_particles_) return;

  glColor3f(0.1f, 0.3f, 0.9f);  // particle color

  const int segments = 16;
  const float r = 0.4f * dx_;  // visual radius

  for (const Particle& p : particles_) {
    float cx = p.x_(0);
    float cy = p.x_(1);

    glBegin(GL_POLYGON);
    for (int k = 0; k < segments; ++k) {
      float angle = 2.0f * M_PI * static_cast<float>(k) / static_cast<float>(segments);
      float px = cx + r * cosf(angle);
      float py = cy + r * sinf(angle);
      glVertex2f(px, py);
    }
    glEnd();
  }
}

FluidSim::Boundary::Boundary(const Vector2s& center, const Vector2s& parameter, 
                            BOUNDARY_TYPE type, bool inside_flag) {
  center_ = center;
  parameter_ = parameter;
  type_ = type;
  sign_ = inside_flag ? -1.0f : 1.0f;
}

Particle::Particle(const Vector2s& x, const Vector2s& v, const scalar& radius, 
                  const scalar& density)
    : x_(x), v_(v), radii_(radius), mass_(static_cast<scalar>((4.0 / 3.0) * M_PI * radius * radius * radius * density)) {
  c_.setZero();     // Affine matrix
  buf0_.setZero();  // Temp position buffer
}

Particle::Particle(const Particle& other) : x_(other.x_), v_(other.v_), radii_(other.radii_), mass_(other.mass_) {
  c_.setZero();
  buf0_.setZero();
}

// Perform iterative extrapolation of valid velocity values using a simple neighborhood averaging method

void extrapolate(Array2s& grid, Array2s& old_grid, const Array2s& grid_weight, const Array2s& grid_liquid_weight, Array2c& valid, Array2c old_valid_,
                 const Vector2i& offset, int num_layers) {
  if (num_layers <= 0) return;

  // Clear the boundary of valid flags
  for (int j = 0; j < valid.nj; ++j) {
    valid(0, j) = 0;
    valid(valid.ni - 1, j) = 0;
  }
  for (int i = 0; i < valid.ni; ++i) {
    valid(i, 0) = 0;
    valid(i, valid.nj - 1) = 0;
  }

  // Prepare alternating buffers for grid and valid flags
  Array2s* grids[2] = {&grid, &old_grid};
  Array2c* valids[2] = {&valid, &old_valid_};

  // Mark initial valid cells inside the domain
  for (int j = 1; j < grid.nj - 1; ++j) {
    for (int i = 1; i < grid.ni - 1; ++i) {
      bool in_fluid = (grid_weight(i, j) > 0);
      bool near_interface = (grid_liquid_weight(i, j) < 0) || (grid_liquid_weight(i + offset(0), j + offset(1)) < 0);
      valid(i, j) = static_cast<char>(in_fluid && near_interface);
    }
  }

  old_valid_ = valid;
  old_grid = grid;

  // Extrapolate in multiple layers
  for (int layer = 0; layer < num_layers; ++layer) {
    Array2s* src_grid = grids[layer % 2];
    Array2s* dst_grid = grids[(layer + 1) % 2];
    Array2c* src_valid = valids[layer % 2];
    Array2c* dst_valid = valids[(layer + 1) % 2];

    parallel_for(1, grid.nj - 1, [&](int j) {
      for (int i = 1; i < grid.ni - 1; ++i) {
        if (!(*src_valid)(i, j)) {
          scalar total = 0.0;
          int contributors = 0;

          // Check 4 neighbors and average valid values
          if ((*src_valid)(i + 1, j)) {
            total += (*src_grid)(i + 1, j);
            ++contributors;
          }
          if ((*src_valid)(i - 1, j)) {
            total += (*src_grid)(i - 1, j);
            ++contributors;
          }
          if ((*src_valid)(i, j + 1)) {
            total += (*src_grid)(i, j + 1);
            ++contributors;
          }
          if ((*src_valid)(i, j - 1)) {
            total += (*src_grid)(i, j - 1);
            ++contributors;
          }

          if (contributors > 0) {
            (*dst_grid)(i, j) = total / static_cast<scalar>(contributors);
            (*dst_valid)(i, j) = 1;
          }
        }
      }
    });

    // Sync data back if more iterations remain
    if (layer + 1 < num_layers) {
      for (int j = 1; j < grid.nj - 1; ++j) {
        for (int i = 1; i < grid.ni - 1; ++i) {
          if (!(*src_valid)(i, j)) {
            (*src_grid)(i, j) = (*dst_grid)(i, j);
            (*src_valid)(i, j) = (*dst_valid)(i, j);
          }
        }
      }
    }
  }
}