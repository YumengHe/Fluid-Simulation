#pragma once
#include <vector>
#include <iostream>

class Fluid_Grid{
public:
    // grid dimension
    int g_width; // number of cells in horizontal direction
    int g_height; // number of cells in vertical direction

    // physics parameter
    int g_num_iteration; // number of iteration for linear solver
    float g_dt; // time step
    float g_diffusion; // coefficient for diffusion
    float g_viscosity; // coefficient for viscosity
    
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
    void initialization(int N, int num_iteration, int dt, float diffusion, float viscosity, std::vector< std::vector<float> >  &velocity_x, std::vector< std::vector<float> >  velocity_y, std::vector< std::vector<float> >  pressure, std::vector< std::vector<float> >  density);
    // the main simulation function
    void simulation();
};

// helper functions
void set_bnd(int N, int b, std::vector< std::vector<float> > &x);
void lin_solve(int b, std::vector< std::vector<float> > &x, std::vector< std::vector<float> > &x0, float a, float c, int num_iteration, int N);
void add_source(int N, std::vector< std::vector<float> > &x, std::vector< std::vector<float> > &s, float dt);

// physics functions
void advect(int N, int b, std::vector< std::vector<float> > d, std::vector< std::vector<float> > d0, std::vector< std::vector<float> > v_x, std::vector< std::vector<float> > v_y, float dt);
void diffuse(int N, int b, std::vector< std::vector<float> > &x, std::vector< std::vector<float> > &x0, float diff, float dt, int num_iteration);
void project(int N, std::vector< std::vector<float> > &velocity_x, std::vector< std::vector<float> > &velocity_y, std::vector< std::vector<float> > &p, std::vector< std::vector<float> > &div, int num_iteration);

// core functions
void dens_step();
void vel_step();

// overload cout to print every element of velocity
std::ostream& operator<<(std::ostream& os, const Fluid_Grid &grid);
