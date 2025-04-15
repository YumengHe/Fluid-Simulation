#pragma once
#include <vector>

class Fluid_Grid{
public:
    // grid dimension
    int g_width; // number of cells in horizontal direction
    int g_height; // number of cells in vertical direction

    // physics parameter
    float g_dt; // time step
    float g_diffusion; // coefficient for diffusion
    float g_viscosity; // coefficient for viscosity
    int g_num_iteration; // number of iteration for linear solver

    //stored at the boundary of the cell
    std::vector<std::vector<float>> g_velocity_x;// Velocity in the x direction
    std::vector<std::vector<float>> g_velocity_y;// Velocity in the y direction

    //stored at the center of the cell
    std::vector<std::vector<float>> g_pressure;
    std::vector<std::vector<float>> g_density;

    //Constructor
    Fluid_Grid(int width,int height, float dt, float diffusion, float viscosity, int num_iteration) {
        
    };

    //Initialize
    void initializeGrid();

    //add core function
    //Step 1: Advection（对密度和速度进行推进）
    void advect();
    //Step 2: Apply Forces
    void diffuse();
    //Step 3: Projection（解 Poisson 方程）  
    void project();
    // apply boundary conditions (e.g. walls)
    void setBoundary();

    void simulation();
};