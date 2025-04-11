#pragma once
#include <vector>

class Fluid_Grid{
public:
    int grid_w; //width of grid (x direction)
    int grid_h; // height of grid (y direction)
    float cell_size;

    //stored at the boundary of the cell
    std::vector<std::vector<float>> velocity_x;// Velocity in the x direction
    std::vector<std::vector<float>> velocity_y;// Velocity in the y direction

    //stored at the center of the cell
    std::vector<std::vector<float>> pressure;
    std::vector<std::vector<float>> density;

    //Constructor
    Fluid_Grid(int width,int height,float cell_size);

    //Initialize
    void initializeGrid();

    //add core function
    //Step 1: Advection（对密度和速度进行推进）
    //Step 2: Apply Forces
    //Step 3: Projection（解 Poisson 方程）  

};