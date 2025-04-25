#pragma once
#include <vector>
#include <iostream>
#include "../../include/Eigen/Dense"

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

    // For each particle, stores the gradients of interpolation weights with respect to x and y directions
    // Used in G2P
    std::vector<std::vector<Eigen::Vector2f>> g_weights_x; // Gradients of interpolation weights with respect to x
    std::vector<std::vector<Eigen::Vector2f>> g_weights_y;// Gradients of interpolation weights with respect to y
     Grid_PIC(int g_width,int g_height,float dx){
        grid_width = g_width;
        grid_height = g_height;
        grid_dx = dx;
    
        g_velocity_x=std::vector<std::vector<float>>(grid_width+1,std::vector<float>(grid_height,0.0f));
        g_velocity_y=std::vector<std::vector<float>>(grid_width,std::vector<float>(grid_height+1,0.0f));

    
        g_pressure=std::vector<std::vector<float>>(grid_width,std::vector<float>(grid_height,0.0f));

        g_mass_x=std::vector<std::vector<float>>(grid_width+1,std::vector<float>(grid_height,0.0f));
        g_mass_y=std::vector<std::vector<float>>(grid_width,std::vector<float>(grid_height+1,0.0f));
                // Initialize APIC weights
        g_weights_x = std::vector<std::vector<Eigen::Vector2f>>(
            grid_width + 1,
            std::vector<Eigen::Vector2f>(grid_height, Eigen::Vector2f::Zero())
        );
        g_weights_y = std::vector<std::vector<Eigen::Vector2f>>(
            grid_width,
            std::vector<Eigen::Vector2f>(grid_height + 1, Eigen::Vector2f::Zero())
        );
    }
};