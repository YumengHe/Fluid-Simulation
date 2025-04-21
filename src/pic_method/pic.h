#pragma once
#include "grid_pic.h"
#include "particle_pic.h"

//initialize particles
void init_particle(std::vector<Particle_PIC>& particles,int num_w,int num_h,float dx);
//Simulation step
void simul_step(std::vector<Particle_PIC>& particles,Grid_PIC& grid,float dt);
