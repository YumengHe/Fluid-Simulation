#include "pic.h"
#include "../../include/Eigen/Sparse"
#include "../../include/Eigen/IterativeLinearSolvers"
using namespace Eigen;

// helper function
// Compute B-spline weights and their gradients for APIC interpolation
// x_f: Fractional part of x coordinate
// y_f: Fractional part of y coordinate
// weights_x: Array to store x-direction weights
// weights_y: Array to store y-direction weights
// dweights_x: Array to store x-direction weight gradients
// dweights_y: Array to store y-direction weight gradients
// dx: Grid cell size
void compute_weights_and_gradients(float x_f, float y_f, 
                                 float* weights_x, float* weights_y,
                                 float* dweights_x, float* dweights_y,
                                 float dx) {
    // B-spline weights
    weights_x[0] = 0.5f * (1.5f - x_f) * (1.5f - x_f);
    weights_x[1] = 0.75f - (x_f - 1.0f) * (x_f - 1.0f);
    weights_x[2] = 0.5f * (x_f - 0.5f) * (x_f - 0.5f);

    weights_y[0] = 0.5f * (1.5f - y_f) * (1.5f - y_f);
    weights_y[1] = 0.75f - (y_f - 1.0f) * (y_f - 1.0f);
    weights_y[2] = 0.5f * (y_f - 0.5f) * (y_f - 0.5f);

    // weight gradients
    dweights_x[0] = (x_f - 1.5f) / dx;
    dweights_x[1] = (-2.0f * (x_f - 1.0f)) / dx;
    dweights_x[2] = (x_f - 0.5f) / dx;

    dweights_y[0] = (y_f - 1.5f) / dx;
    dweights_y[1] = (-2.0f * (y_f - 1.0f)) / dx;
    dweights_y[2] = (y_f - 0.5f) / dx;
}

// APIC formula：p2g
// m_i * v_i = ∑_p (w_ip * m_p * (v_p + C_p * (x_i - x_p)))

// helper function
// Transfer a single particle's velocity component to nearby grid cells using APIC method
// g_velocity: 2D grid storing velocity
// g_mass: 2D grid storing mass weights
// g_weights: 2D grid storing weight gradients for APIC
// p: Particle containing position and affine matrix
// velocity_component: Particle velocity in either x or y
// offset_x: Grid offset in x direction
// offset_y: Grid offset in y direction
// grid_width: Width of the velocity grid
// grid_height: Height of the velocity grid
// dx: Grid cell size
void transfer_velocity_apic(std::vector<std::vector<float>>& g_velocity,
                          std::vector<std::vector<float>>& g_mass,
                          std::vector<std::vector<Eigen::Vector2f>>& g_weights,
                          const Particle_PIC& p,
                          float velocity_component,
                          float offset_x, float offset_y,
                          int grid_width, int grid_height,
                          float dx) {
    float g_x = p.x/dx - offset_x;
    float g_y = p.y/dx - offset_y;

    int int_x = static_cast<int>(g_x);
    int int_y = static_cast<int>(g_y);

    float x_f = g_x - int_x;
    float y_f = g_y - int_y;

    float weights_x[3], weights_y[3];
    float dweights_x[3], dweights_y[3];
    compute_weights_and_gradients(x_f, y_f, weights_x, weights_y, dweights_x, dweights_y, dx);

    // Compute D_p matrix (second moment of weights)
    Matrix2f D_p = Matrix2f::Zero();
    for (int i = 0; i < 3; i++) {
        int x_i = int_x + i;
        if (x_i < 0 || x_i >= grid_width) continue;

        for (int j = 0; j < 3; j++) {
            int y_i = int_y + j;
            if (y_i < 0 || y_i >= grid_height) continue;

            float dx_i = (x_i + offset_x) * dx - p.x;
            float dy_i = (y_i + offset_y) * dx - p.y;
            float weight = weights_x[i] * weights_y[j];
            
            // Accumulate second moments
            D_p(0,0) += dx_i * dx_i * weight;
            D_p(0,1) += dx_i * dy_i * weight;
            D_p(1,0) += dy_i * dx_i * weight;
            D_p(1,1) += dy_i * dy_i * weight;
        }
    }

    // Calculate inverse of D_p matrix
    Matrix2f D_p_inv;
    float det = D_p(0,0) * D_p(1,1) - D_p(0,1) * D_p(1,0);
    if (std::abs(det) > 1e-10f) {
        D_p_inv(0,0) = D_p(1,1) / det;
        D_p_inv(0,1) = -D_p(0,1) / det;
        D_p_inv(1,0) = -D_p(1,0) / det;
        D_p_inv(1,1) = D_p(0,0) / det;
    } else {
        D_p_inv = Matrix2f::Identity();
    }

    // Calculate C_p = B_p * D_p^{-1}
    Matrix2f C_p = p.B * D_p_inv;

    
    // 3x3 grid
    for (int i = 0; i < 3; i++) {
        int x_i = int_x + i;
        if (x_i < 0 || x_i >= grid_width) continue;

        for (int j = 0; j < 3; j++) {
            int y_i = int_y + j;
            if (y_i < 0 || y_i >= grid_height) continue;

            float weight = weights_x[i] * weights_y[j];
            Vector2f dweight(dweights_x[i] * weights_y[j], weights_x[i] * dweights_y[j]);

            // calculate position difference between particle and grid cell
            float dx_i = (x_i + offset_x) * dx - p.x;
            float dy_i = (y_i + offset_y) * dx - p.y;
            
            // affine term: C_p * (x_i - x_p)
            float affine_term = C_p(0,0) * dx_i + C_p(0,1) * dy_i;  // 对于x方向速度
            if (velocity_component == p.velocity_y) {
                affine_term = C_p(1,0) * dx_i + C_p(1,1) * dy_i;  // 对于y方向速度
            }

            // Complete APIC transfer formula: v_p + C_p * (x_i - x_p)
            g_velocity[x_i][y_i] += (velocity_component + affine_term) * weight;
            g_mass[x_i][y_i] += weight;
            g_weights[x_i][y_i] = dweight;
        }
    }
}

// Transfer particle velocities and affine matrices to the grid using APIC method
// particles: Array containing PIC/APIC particles
// grid: MAC grid storing velocities, masses, and APIC specific weights
void p2g_apic(const std::vector<Particle_PIC>& particles, Grid_PIC& grid) {
    reset_grid(grid);
    
    for (const auto& p : particles) {
        // transfer x direction velocity
        transfer_velocity_apic(grid.g_velocity_x, grid.g_mass_x, grid.g_weights_x,
                             p, p.velocity_x, 0.5f, 0.0f,
                             grid.grid_width + 1, grid.grid_height, grid.grid_dx);
        
        // transfer y direction velocity
        transfer_velocity_apic(grid.g_velocity_y, grid.g_mass_y, grid.g_weights_y,
                             p, p.velocity_y, 0.0f, 0.5f,
                             grid.grid_width, grid.grid_height + 1, grid.grid_dx);
    }

    // normalize
    for (int i = 0; i <= grid.grid_width; i++) {
        for (int j = 0; j < grid.grid_height; j++) {
            if (grid.g_mass_x[i][j] > 0.0f) {
                grid.g_velocity_x[i][j] /= grid.g_mass_x[i][j];
            }
        }
    }

    for (int i = 0; i < grid.grid_width; i++) {
        for (int j = 0; j <= grid.grid_height; j++) {
            if (grid.g_mass_y[i][j] > 0.0f) {
                grid.g_velocity_y[i][j] /= grid.g_mass_y[i][j];
            }
        }
    }
}

// Update particle velocities and affine matrices by interpolating from the grid using APIC
// particles: Array of particles to update
// grid: MAC grid containing velocities and APIC specific weights
void g2p_apic(std::vector<Particle_PIC>& particles, const Grid_PIC& grid) {
    for (auto& p : particles) {
        float p_x = p.x;
        float p_y = p.y;
        
        // Calculate interpolation weights and gradients
        float g_x = p_x/grid.grid_dx - 0.5f;
        float g_y = p_y/grid.grid_dx;
        
        int int_x = static_cast<int>(g_x);
        int int_y = static_cast<int>(g_y);
        
        float x_f = g_x - int_x;
        float y_f = g_y - int_y;
        
        float weights_x[3], weights_y[3];
        float dweights_x[3], dweights_y[3];
        compute_weights_and_gradients(x_f, y_f, weights_x, weights_y, dweights_x, dweights_y, grid.grid_dx);
        
        // Update particle velocity and affine matrix
        p.velocity_x = 0.0f;
        p.velocity_y = 0.0f;
        p.B.setZero();
        
        // Calculate x-direction velocity and affine matrix
        for (int i = 0; i < 3; i++) {
            int x_i = int_x + i;
            if (x_i < 0 || x_i >= grid.grid_width + 1) continue;
            
            for (int j = 0; j < 3; j++) {
                int y_i = int_y + j;
                if (y_i < 0 || y_i >= grid.grid_height) continue;
                
                float weight = weights_x[i] * weights_y[j];
                Vector2f grid_pos = Vector2f((x_i + 0.5f) * grid.grid_dx, y_i * grid.grid_dx);
                Vector2f delta = grid_pos - Vector2f(p.x, p.y);
                Vector2f v_grid(grid.g_velocity_x[x_i][y_i], 0.0f);
                
                p.velocity_x += grid.g_velocity_x[x_i][y_i] * weight;
                p.B += v_grid * delta.transpose() * weight;
            }
        }
        
        // Calculate y-direction velocity and affine matrix
        g_x = p_x/grid.grid_dx;
        g_y = p_y/grid.grid_dx - 0.5f;
        int_x = static_cast<int>(g_x);
        int_y = static_cast<int>(g_y);
        x_f = g_x - int_x;
        y_f = g_y - int_y;
        
        compute_weights_and_gradients(x_f, y_f, weights_x, weights_y, dweights_x, dweights_y, grid.grid_dx);
        
        for (int i = 0; i < 3; i++) {
            int x_i = int_x + i;
            if (x_i < 0 || x_i >= grid.grid_width) continue;
            
            for (int j = 0; j < 3; j++) {
                int y_i = int_y + j;
                if (y_i < 0 || y_i >= grid.grid_height + 1) continue;
                
                float weight = weights_x[i] * weights_y[j];
                Vector2f grid_pos = Vector2f(x_i * grid.grid_dx, (y_i + 0.5f) * grid.grid_dx);
                Vector2f delta = grid_pos - Vector2f(p.x, p.y);
                Vector2f v_grid(0.0f, grid.g_velocity_y[x_i][y_i]);
                
                p.velocity_y += grid.g_velocity_y[x_i][y_i] * weight;
                p.B += v_grid * delta.transpose() * weight;
            }
        }
    }
}

// Perform a complete APIC simulation step
// particles: Array of particles in the simulation
// grid: MAC grid storing velocities, pressure, and APIC specific data
// dt: Time step size
void simul_step_apic(std::vector<Particle_PIC>& particles, Grid_PIC& grid, float dt) {
    p2g_apic(particles, grid);
    apply_gravity(grid, dt);
    solve_pressure(grid, dt);
    g2p_apic(particles, grid);
    advect(particles, dt, grid);
}