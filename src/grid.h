#pragma once
#include <vector>
#include <iostream>

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
    std::vector< std::vector<float> > g_velocity_x;// Velocity in the x direction
    std::vector< std::vector<float> > g_velocity_y;// Velocity in the y direction
    std::vector< std::vector<float> > g_velocity_x0;// Old velocity in the x direction
    std::vector< std::vector<float> > g_velocity_y0;// Old velocity in the y direction

    //stored at the center of the cell
    std::vector< std::vector<float> > g_pressure;
    std::vector< std::vector<float> > g_density;
    std::vector< std::vector<float> > g_density0;

    //Constructor
    Fluid_Grid();
    Fluid_Grid(int width,int height, float dt, float diffusion, float viscosity, int num_iteration);

    // initialize grid
    void initialization();
    // the main simulation function
    void simulation();

    // core functions
    // Step 1: Advection（对密度和速度进行推进
    // Step 2: Apply Forces
    // Step 3: Projection（解 Poisson 方程）
};

// overload cout to print every element of velocity
std::ostream& operator<<(std::ostream& os, const Fluid_Grid &grid) {
    os << "g_velocity_x:\n";
    for (int j = 0; j < grid.g_height; j++) {
        for (int i = 0; i < grid.g_width; i++) {
            os << grid.g_velocity_x[j][i] << " ";
        }
        os << "\n";
    }

    os << "\n";

    os << "g_velocity_y:\n";
    for (int j = 0; j < grid.g_height; j++) {
        for (int i = 0; i < grid.g_width; i++) {
            os << grid.g_velocity_y[j][i] << " ";
        }
        os << "\n";
    }

    return os;
}