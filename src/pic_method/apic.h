#pragma once
#include "grid_pic.h"
#include "particle_pic.h"
#include "../../include/Eigen/Dense"

// Compute B-spline weights and gradients
void compute_weights_and_gradients(float x_f, float y_f, 
                                 float* weights_x, float* weights_y,
                                 float* dweights_x, float* dweights_y,
                                 float dx);

// Transfer particle velocity to grid using APIC method
void transfer_velocity_apic(std::vector<std::vector<float>>& g_velocity,
                          std::vector<std::vector<float>>& g_mass,
                          std::vector<std::vector<Eigen::Vector2f>>& g_weights,
                          const Particle_PIC& p,
                          float velocity_component,
                          float offset_x, float offset_y,
                          int grid_width, int grid_height,
                          float dx);

// Complete particle-to-grid transfer process
void p2g_apic(const std::vector<Particle_PIC>& particles, Grid_PIC& grid);

// Grid-to-particle transfer process
void g2p_apic(std::vector<Particle_PIC>& particles, const Grid_PIC& grid);

// Complete APIC simulation step
void simul_step_apic(std::vector<Particle_PIC>& particles, Grid_PIC& grid, float dt);

extern std::vector<Particle_PIC> apic_particles;
extern Grid_PIC apic_grid;
extern float dt;