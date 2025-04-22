#pragma once
#include <vector>
#include <iostream>

//MAC grid structure
struct Grid_PIC{
    int grid_width;//number of grids in x direction
    int grid_height;//number of grids in y direction
    float grid_dx;//size of single grid

    std::vector<std::vector<float>> g_velocity_x;//velocity in the x direction,Stored on the left boundary of the cell
    std::vector<std::vector<float>> g_velocity_y;//velocity in the y direction,Stored on the bottom boundary of the cell

    std::vector<std::vector<float>> g_pressure;

     // mass weights from P2G
     std::vector<std::vector<float>> g_mass_x;//mass weight in the x direction
     std::vector<std::vector<float>> g_mass_y;//mass weight in the y direction

     Grid_PIC(int g_width,int g_height,float dx){
        grid_width = g_width;
        grid_height = grid_height;
        grid_dx = dx;
    
        g_velocity_x=std::vector<std::vector<float>>(grid_width+1,std::vector<float>(grid_height,0.0f));
        g_velocity_y=std::vector<std::vector<float>>(grid_width,std::vector<float>(grid_height+1,0.0f));

    
        g_pressure=std::vector<std::vector<float>>(grid_width,std::vector<float>(grid_height,0.0f));

        g_mass_x=std::vector<std::vector<float>>(grid_width+1,std::vector<float>(grid_height,0.0f));
        g_mass_y=std::vector<std::vector<float>>(grid_width,std::vector<float>(grid_height+1,0.0f));
    }
};