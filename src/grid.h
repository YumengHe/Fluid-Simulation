#pragma once
#include <vector>
#include <include/eigen-3.4.0/Eigen/Dense>

class Fluid_Grid{
public:
    // grid dimension
    int grid_w; // number of cells in horizontal direction
    int grid_h; // number of cells in vertical direction

    // physics parameter
    float dt; // time step
    float diffusion; // coefficient for diffusion
    float viscosity; // coefficient for viscosity
    int num_iteration; // number of iteration for linear solver

    //stored at the boundary of the cell
    std::vector<std::vector<float>> velocity_x;// Velocity in the x direction
    std::vector<std::vector<float>> velocity_y;// Velocity in the y direction

    //stored at the center of the cell
    std::vector<std::vector<float>> pressure;
    std::vector<std::vector<float>> density;

    //Constructor
    Fluid_Grid(int width,int height) {
        
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