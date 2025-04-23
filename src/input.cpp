#include <GLUT/glut.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "input.h"
#include "grid.h"
#include "particle.h"  // 首先需要包含头文件
extern std::vector<Particle> particles;  // 声明这是一个外部变量

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
}

void idle()
{
    if (current_grid) {
        current_grid->simulation();  // run simulation for grid-based fluid
    }
    else if (!particles.empty()) {
        Update();  // run simulation for particle-based fluid
    }
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

    // Read density field
    grid.g_density.resize(height);
    for (int i = 0; i < height; ++i) {
        grid.g_density[i].resize(width);
        for (int j = 0; j < width; ++j) {
            file >> grid.g_density[i][j];
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
        particles.push_back(Particle(x, y));
    }
    
    std::cout << "Loaded " << particles.size() << " particles from " << filename << std::endl;
}