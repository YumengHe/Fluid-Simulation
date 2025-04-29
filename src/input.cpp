#include <GLUT/glut.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "input.h"
#include "grid.h"
#include "particle.h"  // 首先需要包含头文件
#include "./pic_method/particle_pic.h"
#include "./pic_method/grid_pic.h"
#include "constants.h"
#include "./pic_method/apic.h"
#include "./pic_method/pic.h"

int pause = 1;

void mouseClick(int button, int state, int x, int y)
{
    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
    {
        // convert screen coord (x, y) to OpenGL 0~1 region
        float xf = (float)x / glutGet(GLUT_WINDOW_WIDTH);
        float yf = 1.0f - (float)y / glutGet(GLUT_WINDOW_HEIGHT); // flip y coord

        printf("🖱️ Click at: screen = (%d, %d), normalized = (%.3f, %.3f)\n", x, y, xf, yf);
    }
}

void handleKeypress(unsigned char key, int x, int y)
{
    if (key == 27){            // 27 is the ASCII value of ESC
        exit(0); // esc
    }

    if (key == 'p') {
        pause = 1 - pause;
    }
}

int frame_counter = 0;
void idle() {
    if (pause == 0) {
        if (current_grid) {
            // std::cout << "doing grid simulation" << std::endl;
            current_grid->simulation();  // run simulation for grid-based fluid
        }
        else if (!particles.empty()) {
            Update();  // run simulation for particle-based fluid
        }else if (!apic_particles.empty()) {
            simul_step_apic(apic_particles, apic_grid, dt);
        }
        else if (!pic_particles.empty()) {
            simul_step(pic_particles, pic_grid, dt_pic);  // 添加PIC的模拟步骤
        }
    }
    

    // for (size_t i = 0; i < apic_particles.size(); ++i) {
    //     std::cout << "Particle " << i 
    //                 << " : (x=" << apic_particles[i].x 
    //                 << ", y=" << apic_particles[i].y << ")\n";
    // }

    glutPostRedisplay(); // request a redraw when cpu is idle
}

bool endsWith(const std::string &str, const std::string &suffix)
{
    return str.size() >= suffix.size() &&
           str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

void loadGrid(const std::string &filename, Fluid_Grid &grid){
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Failed to open grid file: " << filename << std::endl;
        return;
    }

    // Read basic parameters
    int width, height, num_iteration;
    float dt, diffusion, viscosity;
    
    file >> width >> height >> num_iteration;
    file >> dt >> diffusion >> viscosity;
    
    // Initialize grid dimensions and parameters
    grid.g_width = width;
    grid.g_height = height;
    grid.g_num_iteration = num_iteration;
    grid.g_dt = dt;
    grid.g_diffusion = diffusion;
    grid.g_viscosity = viscosity;

    // Initialize velocity and pressure to 0
    grid.g_velocity_x.assign(height, std::vector<float>(width, 0.0f));
    grid.g_velocity_y.assign(height, std::vector<float>(width, 0.0f));
    grid.g_pressure.assign(height, std::vector<float>(width, 0.0f));

    /*
    // Read velocity field x component
    grid.g_velocity_x.resize(height);
    for (int i = 0; i < height; ++i) {
        grid.g_velocity_x[i].resize(width);
        for (int j = 0; j < width; ++j) {
            file >> grid.g_velocity_x[i][j];
        }
    }

    // Read velocity field y component
    grid.g_velocity_y.resize(height);
    for (int i = 0; i < height; ++i) {
        grid.g_velocity_y[i].resize(width);
        for (int j = 0; j < width; ++j) {
            file >> grid.g_velocity_y[i][j];
        }
    }

    // Read pressure field
    grid.g_pressure.resize(height);
    for (int i = 0; i < height; ++i) {
        grid.g_pressure[i].resize(width);
        for (int j = 0; j < width; ++j) {
            file >> grid.g_pressure[i][j];
        }
    }
    */

    // Read density field
    grid.g_density.resize(height);
    for (int i = 0; i < height; ++i) {
        grid.g_density[i].resize(width);
        for (int j = 0; j < width; ++j) {
            file >> grid.g_density[i][j];
            if (grid.g_density[i][j] != 0) {
                std::cout << "Cell [" << i << "][" << j << "] has density :" << 
                grid.g_density[i][j] << std::endl;
            }
        }
    }
    grid.initialization(width, num_iteration, dt, diffusion, viscosity,
                        grid.g_velocity_x, grid.g_velocity_y,
                        grid.g_pressure, grid.g_density);
    file.close();
    
    std::cout << "Loaded grid configuration from: " << filename 
              << "\nGrid size: " << width << "x" << height 
              << "\nIterations: " << num_iteration << std::endl;
}

void loadParticles(const std::string& filename) {
    particles.clear();
    std::ifstream file(filename);
    std::string line;
    float x, y;
    while (file >> x >> y) {
        float jitter = static_cast<float>(arc4random()) / static_cast<float>(RAND_MAX);
        particles.push_back(Particle(x + jitter, y));
    }
    
    std::cout << "Loaded " << particles.size() << " particles from " << filename << std::endl;
}

void loadAPIC(const std::string &filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }

    // Read time step
    file >> dt;

    // Read grid information
    int grid_width_num, grid_height_num;
    file >> grid_width_num >> grid_height_num;
    float grid_size =  1.0f / grid_width_num;
    apic_grid = Grid_PIC(grid_width_num, grid_height_num, grid_size);
    // Read all particle data
    float x, y, vx, vy;
    float B11, B12, B21, B22;
    while (file >> x >> y >> vx >> vy >> B11 >> B12 >> B21 >> B22) {
        Particle_PIC particle(x / VIEW_WIDTH, y / VIEW_WIDTH, vx, vy);
        particle.B(0,0) = B11;
        particle.B(0,1) = B12;
        particle.B(1,0) = B21;
        particle.B(1,1) = B22;
        // B matrix is already initialized to zero in constructor
        apic_particles.push_back(particle);
    }

    file.close();

    std::cout << "Loaded " << apic_particles.size() << " apic_particles\n";
    std::cout << "Grid size: " << grid_width_num << "x" << grid_height_num << "\n";
    std::cout << "Grid spacing: " << grid_size << "\n";
    std::cout << "dt = " << dt << "\n"; // ✅ 打印一下检查
    std::cout << "First particle: (" << apic_particles[0].x << ", " << apic_particles[0].y << "), "
          << "vx: " << apic_particles[0].velocity_x << ", vy: " << apic_particles[0].velocity_y << std::endl;

}

void loadPIC(const std::string &filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }

    // 读取时间步长
    file >> dt_pic;

    // 读取网格信息
    int grid_width_num, grid_height_num;
    file >> grid_width_num >> grid_height_num;
    float grid_size = 1 / grid_width_num;
    pic_grid = Grid_PIC(grid_width_num, grid_height_num, grid_size);

    // 清空现有粒子数据
    pic_particles.clear();

    // 读取粒子数据
    float x, y, vx, vy;
    while (file >> x >> y >> vx >> vy) {
        Particle_PIC particle(x / VIEW_WIDTH, y / VIEW_WIDTH, vx, vy);
        pic_particles.push_back(particle);
    }

    file.close();

    std::cout << "Loaded " << pic_particles.size() << " PIC particles\n";
    std::cout << "Grid size: " << grid_width_num << "x" << grid_height_num << "\n";
    std::cout << "Grid spacing: " << grid_size << "\n";
    std::cout << "dt_pic = " << dt_pic << "\n";
}