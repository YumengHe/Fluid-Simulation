#include "pic.h"
#include "../../include/Eigen/Sparse"
#include "../../include/Eigen/IterativeLinearSolvers"
using namespace Eigen;
#include <iomanip> // for std::setprecision


// Define global variables
std::vector<Particle_PIC> apic_particles;
Grid_PIC apic_grid(0, 0, 0.0f);  // 初始化为默认值
float dt = 0.01f;  // 设置默认时间步长

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

    // Calculate inverse of D_p matrix with safety
    Matrix2f D_p_inv;
    float det = D_p(0,0) * D_p(1,1) - D_p(0,1) * D_p(1,0);
    if (std::abs(det) > 1e-6f) { // 比原来更宽松，防止小数值炸掉
        D_p_inv(0,0) = D_p(1,1) / det;
        D_p_inv(0,1) = -D_p(0,1) / det;
        D_p_inv(1,0) = -D_p(1,0) / det;
        D_p_inv(1,1) = D_p(0,0) / det;
    } else {
        D_p_inv = Matrix2f::Identity(); // 保护
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
            float affine_term = (velocity_component == p.velocity_y) ?
                (C_p(1,0) * dx_i + C_p(1,1) * dy_i) :
                (C_p(0,0) * dx_i + C_p(0,1) * dy_i);

            float v_transfer = velocity_component + affine_term;

            // Check for NaN
            if (std::isnan(v_transfer) || std::isnan(weight)) continue;

            // Complete APIC transfer formula: v_p + C_p * (x_i - x_p)
            g_velocity[x_i][y_i] += v_transfer * weight;
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
    float dx = grid.grid_dx;

    for (const auto& p : particles) {
        // For velocity_x (offset 0.5,0)
        {
            float g_x = p.x / dx - 0.5f;
            float g_y = p.y / dx;
            int base_x = (int)g_x;
            int base_y = (int)g_y;
            float fx = g_x - base_x;
            float fy = g_y - base_y;

            float weights_x[3], weights_y[3];
            float dweights_x[3], dweights_y[3];
            compute_weights_and_gradients(fx, fy, weights_x, weights_y, dweights_x, dweights_y, dx);

            for (int i = 0; i < 3; ++i) {
                int ix = base_x + i;
                if (ix < 0 || ix >= grid.grid_width + 1) continue;
                for (int j = 0; j < 3; ++j) {
                    int iy = base_y + j;
                    if (iy < 0 || iy >= grid.grid_height) continue;

                    float weight = weights_x[i] * weights_y[j];
                    float dx_i = (ix + 0.5f) * dx - p.x;
                    float dy_i = (iy + 0.0f) * dx - p.y;
                    float affine = p.B(0,0) * dx_i + p.B(0,1) * dy_i;

                    grid.g_velocity_x[ix][iy] += (p.velocity_x + affine) * weight;
                    grid.g_mass_x[ix][iy] += weight;
                }
            }
        }

        // For velocity_y (offset 0,0.5)
        {
            float g_x = p.x / dx;
            float g_y = p.y / dx - 0.5f;
            int base_x = (int)g_x;
            int base_y = (int)g_y;
            float fx = g_x - base_x;
            float fy = g_y - base_y;

            float weights_x[3], weights_y[3];
            float dweights_x[3], dweights_y[3];
            compute_weights_and_gradients(fx, fy, weights_x, weights_y, dweights_x, dweights_y, dx);

            for (int i = 0; i < 3; ++i) {
                int ix = base_x + i;
                if (ix < 0 || ix >= grid.grid_width) continue;
                for (int j = 0; j < 3; ++j) {
                    int iy = base_y + j;
                    if (iy < 0 || iy >= grid.grid_height + 1) continue;

                    float weight = weights_x[i] * weights_y[j];
                    float dx_i = (ix + 0.0f) * dx - p.x;
                    float dy_i = (iy + 0.5f) * dx - p.y;
                    float affine = p.B(1,0) * dx_i + p.B(1,1) * dy_i;

                    grid.g_velocity_y[ix][iy] += (p.velocity_y + affine) * weight;
                    grid.g_mass_y[ix][iy] += weight;
                }
            }
        }
    }

    // Normalize by mass
    for (int i = 0; i <= grid.grid_width; ++i) {
        for (int j = 0; j < grid.grid_height; ++j) {
            if (grid.g_mass_x[i][j] > 1e-8f) {
                grid.g_velocity_x[i][j] /= grid.g_mass_x[i][j];
            }
        }
    }
    for (int i = 0; i < grid.grid_width; ++i) {
        for (int j = 0; j <= grid.grid_height; ++j) {
            if (grid.g_mass_y[i][j] > 1e-8f) {
                grid.g_velocity_y[i][j] /= grid.g_mass_y[i][j];
            }
        }
    }
}

// Update particle velocities and affine matrices by interpolating from the grid using APIC
// particles: Array of particles to update
// grid: MAC grid containing velocities and APIC specific weights
void g2p_apic(std::vector<Particle_PIC>& particles, const Grid_PIC& grid) {
    float dx = grid.grid_dx;

    for (auto& p : particles) {
        p.new_velocity_x = 0.0f;
        p.new_velocity_y = 0.0f;
        p.B.setZero(); // Reset affine matrix

        // For velocity_x (staggered at i+0.5,j)
        {
            float g_x = p.x / dx - 0.5f;
            float g_y = p.y / dx;
            int base_x = (int)g_x;
            int base_y = (int)g_y;
            float fx = g_x - base_x;
            float fy = g_y - base_y;

            float weights_x[3], weights_y[3];
            float dweights_x[3], dweights_y[3];
            compute_weights_and_gradients(fx, fy, weights_x, weights_y, dweights_x, dweights_y, dx);

            for (int i = 0; i < 3; ++i) {
                int ix = base_x + i;
                if (ix < 0 || ix >= grid.grid_width + 1) continue;
                for (int j = 0; j < 3; ++j) {
                    int iy = base_y + j;
                    if (iy < 0 || iy >= grid.grid_height) continue;

                    float weight = weights_x[i] * weights_y[j];
                    float v = grid.g_velocity_x[ix][iy];

                    p.new_velocity_x += v * weight;

                    Eigen::Vector2f dweight(dweights_x[i] * weights_y[j], weights_x[i] * dweights_y[j]);
                    p.B(0,0) += v * dweight.x();
                    p.B(0,1) += v * dweight.y();
                }
            }
        }

        // For velocity_y (staggered at i,j+0.5)
        {
            float g_x = p.x / dx;
            float g_y = p.y / dx - 0.5f;
            int base_x = (int)g_x;
            int base_y = (int)g_y;
            float fx = g_x - base_x;
            float fy = g_y - base_y;

            float weights_x[3], weights_y[3];
            float dweights_x[3], dweights_y[3];
            compute_weights_and_gradients(fx, fy, weights_x, weights_y, dweights_x, dweights_y, dx);

            for (int i = 0; i < 3; ++i) {
                int ix = base_x + i;
                if (ix < 0 || ix >= grid.grid_width) continue;
                for (int j = 0; j < 3; ++j) {
                    int iy = base_y + j;
                    if (iy < 0 || iy >= grid.grid_height + 1) continue;

                    float weight = weights_x[i] * weights_y[j];
                    float v = grid.g_velocity_y[ix][iy];

                    p.new_velocity_y += v * weight;

                    Eigen::Vector2f dweight(dweights_x[i] * weights_y[j], weights_x[i] * dweights_y[j]);
                    p.B(1,0) += v * dweight.x();
                    p.B(1,1) += v * dweight.y();
                }
            }
        }

        // Update particle velocity (new_velocity) ✅
        p.velocity_x = p.new_velocity_x;
        p.velocity_y = p.new_velocity_y;
    }
}

// Perform a complete APIC simulation step
// particles: Array of particles in the simulation
// grid: MAC grid storing velocities, pressure, and APIC specific data
// dt: Time step size
// void simul_step_apic(std::vector<Particle_PIC>& particles, Grid_PIC& grid, float dt) {
//     p2g_apic(particles, grid);
//     apply_gravity(grid, dt);
//     solve_pressure(grid, dt);
//     g2p_apic(particles, grid);
//     advect(particles, dt, grid);
// }


void check_particles(const std::vector<Particle_PIC>& particles, const std::string& stage) {
    bool has_nan = false;
    float min_x = 1e10f, max_x = -1e10f;
    float min_y = 1e10f, max_y = -1e10f;
    float min_vx = 1e10f, max_vx = -1e10f;
    float min_vy = 1e10f, max_vy = -1e10f;

    for (const auto& p : particles) {
        if (std::isnan(p.x) || std::isnan(p.y) || std::isnan(p.velocity_x) || std::isnan(p.velocity_y)) {
            has_nan = true;
            break;
        }
        min_x = std::min(min_x, p.x);
        max_x = std::max(max_x, p.x);
        min_y = std::min(min_y, p.y);
        max_y = std::max(max_y, p.y);
        min_vx = std::min(min_vx, p.velocity_x);
        max_vx = std::max(max_vx, p.velocity_x);
        min_vy = std::min(min_vy, p.velocity_y);
        max_vy = std::max(max_vy, p.velocity_y);
    }

    std::cout << "Stage [" << stage << "] Check:\n";
    if (has_nan) {
        std::cout << "⚠️  NaN detected in particles (position or velocity)!\n";
    } else {
        std::cout << std::fixed << std::setprecision(2);
        std::cout << "Particle X range: [" << min_x << ", " << max_x << "]\n";
        std::cout << "Particle Y range: [" << min_y << ", " << max_y << "]\n";
        std::cout << "Particle Vx range: [" << min_vx << ", " << max_vx << "]\n";
        std::cout << "Particle Vy range: [" << min_vy << ", " << max_vy << "]\n";
    }
    std::cout << "------------------------------------\n";
}

void advect_apic(std::vector<Particle_PIC>& particles, const Grid_PIC& grid, float dt) {
    float dx = grid.grid_dx;
    float min_x = 0.0f;
    float max_x = dx * (grid.grid_width);
    float min_y = 0.0f;
    float max_y = dx * (grid.grid_height);

    for (auto& p : particles) {
        // Step 1: 半步，根据粒子自己的 velocity 做
        float mid_x = p.x + 0.5f * dt * p.velocity_x;
        float mid_y = p.y + 0.5f * dt * p.velocity_y;

        // Step 2: 在 mid 位置插值新的速度
        float vx1 = interp_velocity(grid.g_velocity_x, mid_x, mid_y, 0.5f, 0.0f, dx, grid.grid_width + 1, grid.grid_height);
        float vy1 = interp_velocity(grid.g_velocity_y, mid_x, mid_y, 0.0f, 0.5f, dx, grid.grid_width, grid.grid_height + 1);

        // Step 3: 全步更新位置
        p.x += vx1 * dt;
        p.y += vy1 * dt;

        // Step 4: 更新粒子的 velocity
        p.velocity_x = vx1;
        p.velocity_y = vy1;

        // Step 5: 边界处理
        if (p.x < min_x) { p.x = min_x; p.velocity_x = 0.0f; }
        if (p.x > max_x) { p.x = max_x; p.velocity_x = 0.0f; }
        if (p.y < min_y) { p.y = min_y; p.velocity_y = 0.0f; }
        if (p.y > max_y) { p.y = max_y; p.velocity_y = 0.0f; }
    }
}

void particle_collision_apic(std::vector<Particle_PIC>& particles, const Grid_PIC& grid) {
    float dx = grid.grid_dx;
    float padding = 2.0f * dx;
    float min_x = padding;
    float max_x = dx * (grid.grid_width) - padding;
    float min_y = padding;
    float max_y = dx * (grid.grid_height) - padding;

    for (auto& p : particles) {
        if (p.x < min_x) {
            p.x = min_x;
            p.velocity_x *= -0.5f;
        }
        if (p.x > max_x) {
            p.x = max_x;
            p.velocity_x *= -0.5f;
        }
        if (p.y < min_y) {
            p.y = min_y;
            p.velocity_y *= -0.5f;
        }
        if (p.y > max_y) {
            p.y = max_y;
            p.velocity_y *= -0.5f;
        }
    }
}

void check_grid(const Grid_PIC& grid, const std::string& stage) {
    bool has_nan = false;
    float min_v = 1e10f, max_v = -1e10f;

    // check g_velocity_x
    for (int i = 0; i <= grid.grid_width; ++i) {
        for (int j = 0; j < grid.grid_height; ++j) {
            if (!std::isfinite(grid.g_velocity_x[i][j])) {
                has_nan = true;
            }
            min_v = std::min(min_v, grid.g_velocity_x[i][j]);
            max_v = std::max(max_v, grid.g_velocity_x[i][j]);
        }
    }

    // check g_velocity_y
    for (int i = 0; i < grid.grid_width; ++i) {
        for (int j = 0; j <= grid.grid_height; ++j) {
            if (!std::isfinite(grid.g_velocity_y[i][j])) {
                has_nan = true;
            }
            min_v = std::min(min_v, grid.g_velocity_y[i][j]);
            max_v = std::max(max_v, grid.g_velocity_y[i][j]);
        }
    }

    std::cout << "Stage [" << stage << "] Check Grid:\n";
    if (has_nan) {
        std::cout << "⚠️  NaN detected in grid velocities!\n";
    } else {
        std::cout << std::fixed << std::setprecision(2);
        std::cout << "Grid Velocity range: [" << min_v << ", " << max_v << "]\n";
    }
    std::cout << "------------------------------------\n";
}


void clamp_velocity(Grid_PIC& grid) {
    int w = grid.grid_width;
    int h = grid.grid_height;

    // 左右边界：x方向速度设0
    for (int j = 0; j < h; j++) {
        grid.g_velocity_x[0][j] = 0.0f;
        grid.g_velocity_x[w][j] = 0.0f;
    }

    // 上下边界：y方向速度设0
    for (int i = 0; i < w; i++) {
        grid.g_velocity_y[i][0] = 0.0f;
        grid.g_velocity_y[i][h] = 0.0f;
    }
}

void solve_pressure_apic(Grid_PIC& grid, float dt) {
    using SPM = Eigen::SparseMatrix<float>;
    using Triplet = Eigen::Triplet<float>;
    using VEC = Eigen::VectorXf;

    int w = grid.grid_width;
    int h = grid.grid_height;
    int N = w * h;
    float dx = grid.grid_dx;

    SPM A(N, N);
    std::vector<Triplet> coef;
    VEC B(N);

    // Build matrix A and vector B
    for (int i = 0; i < w; i++) {
        for (int j = 0; j < h; j++) {
            int idx = transfer_index(i, j, h);

            // Solid walls：边界格子处理
            if (i == 0 || i == w-1 || j == 0) { 
                coef.emplace_back(idx, idx, 1.0f);
                B[idx] = 0.0f;
                continue;
            }

            coef.emplace_back(idx, idx, 4.0f);
            coef.emplace_back(idx, transfer_index(i-1, j, h), -1.0f);
            coef.emplace_back(idx, transfer_index(i+1, j, h), -1.0f);
            coef.emplace_back(idx, transfer_index(i, j-1, h), -1.0f);
            coef.emplace_back(idx, transfer_index(i, j+1, h), -1.0f);

            // Compute divergence (∇·u)
            float div = (grid.g_velocity_x[i+1][j] - grid.g_velocity_x[i][j]) / dx +
                        (grid.g_velocity_y[i][j+1] - grid.g_velocity_y[i][j]) / dx;
            B[idx] = div;
        }
    }

    A.setFromTriplets(coef.begin(), coef.end());

    Eigen::ConjugateGradient<SPM, Eigen::Lower|Eigen::Upper> solver;
    solver.compute(A);
    VEC x = solver.solve(B);

    if (solver.info() != Eigen::Success) {
        x.setZero();
    }

    // Fill pressure
    for (int i = 0; i < w; i++) {
        for (int j = 0; j < h; j++) {
            grid.g_pressure[i][j] = x[transfer_index(i, j, h)];
        }
    }

    // Correct velocity
    for (int i = 1; i < w; i++) {
        for (int j = 0; j < h; j++) {
            grid.g_velocity_x[i][j] -= (grid.g_pressure[i][j] - grid.g_pressure[i-1][j]) / dx;
        }
    }
    for (int i = 0; i < w; i++) {
        for (int j = 1; j < h; j++) {
            grid.g_velocity_y[i][j] -= (grid.g_pressure[i][j] - grid.g_pressure[i][j-1]) / dx;
        }
    }
}

void apply_gravity_grid(Grid_PIC& grid, float gravity, float dt) {
    for (int i = 0; i < grid.grid_width; i++) {
        for (int j = 0; j <= grid.grid_height; j++) {
            grid.g_velocity_y[i][j] += gravity * dt;  // y方向有重力
        }
    }
}

void apply_damping(Grid_PIC& grid, float damping) {
    for (int i = 0; i <= grid.grid_width; i++) {
        for (int j = 0; j < grid.grid_height; j++) {
            grid.g_velocity_x[i][j] *= damping;
        }
    }
    for (int i = 0; i < grid.grid_width; i++) {
        for (int j = 0; j <= grid.grid_height; j++) {
            grid.g_velocity_y[i][j] *= damping;
        }
    }
}



void simul_step_apic(std::vector<Particle_PIC>& particles, Grid_PIC& grid, float dt) {
    // reset_grid(grid);

    p2g_apic(particles, grid);

    apply_gravity_grid(grid, -9.8f, dt);
    solve_pressure_apic(grid, dt);
    clamp_velocity(grid);
    apply_damping(grid, 0.99f);
    g2p_apic(particles, grid);   

    advect_apic(particles, grid, dt);   // ✅ 用grid里的速度advect
    particle_collision_apic(particles, grid);
    // g2p_apic(particles, grid);
    // 然后 g2p！ 让粒子拿到新的velocity
}