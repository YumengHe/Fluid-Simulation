#pragma once
#include "grid_pic.h"
#include "particle_pic.h"

//initialize particles
void init_particle(std::vector<Particle_PIC>& particles,int num_w,int num_h,float dx);
//Simulation step
void simul_step(std::vector<Particle_PIC>& particles,Grid_PIC& grid,float dt);
//reset grid
void reset_grid(Grid_PIC& grid);
//particle to grid
void p2g(const std::vector<Particle_PIC>& particles,Grid_PIC& grid);
//apply gravity to grid
void apply_gravity(Grid_PIC& grid,float dt);
//solve pressure, projection
void solve_pressure(Grid_PIC& grid,float dt);
//grid to particle
void g2p(std::vector<Particle_PIC>& particles, const Grid_PIC& grid);
//advection
void advect(std::vector<Particle_PIC>& particles, float dt,const Grid_PIC& grid);
float interp_velocity(const std::vector<std::vector<float>>& g_velocity, 
                      float p_x, float p_y,
                      float offset_x, float offset_y,
                      float dx, int grid_width, int grid_height);
int transfer_index(int i,int j, int height);
extern std::vector<Particle_PIC> pic_particles;
extern Grid_PIC pic_grid;
extern float dt_pic;
