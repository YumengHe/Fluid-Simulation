#include <grid.h>
#include <cmath>

// ------------------------------------------------------------
// Constructor
// ------------------------------------------------------------
Fluid_Grid::Fluid_Grid() {

}

Fluid_Grid::Fluid_Grid(int width,int height, float dt, float diffusion, float viscosity, int num_iteration) {
  g_width = width;
  g_height = height;
  g_dt = dt;
  g_diffusion = diffusion;
  g_viscosity = viscosity;
  g_num_iteration = num_iteration;
}

// ------------------------------------------------------------
// Advection
// ------------------------------------------------------------
void Fluid_Grid::advect() {

}

// ------------------------------------------------------------
// Diffusion
// ------------------------------------------------------------
void Fluid_Grid::diffuse() {
  
}

// ------------------------------------------------------------
// Initialization & Simulation
// ------------------------------------------------------------
void Fluid_Grid::initialization(){
  Fluid_Grid grid;
  grid.g_width = 6 + 2;
  grid.g_height = 5 + 2;
  grid.g_dt = 0.1;
  grid.g_diffusion = 10;
  grid.g_viscosity = 10;
  grid.g_num_iteration = 3;

  // resize
  // g_velocity[row/width][column/height]
  g_velocity_x.resize(g_height, std::vector<float>(g_width, 0.0f));
  g_velocity_y.resize(g_height, std::vector<float>(g_width, 0.0f));
  g_velocity_x0.resize(g_height, std::vector<float>(g_width, 0.0f));
  g_velocity_y0.resize(g_height, std::vector<float>(g_width, 0.0f));

  for (int i = 0; i < g_width; i ++) {
    for (int j = 0; j < g_height; j ++) {
      g_velocity_x[j][i] = 0.0;
      g_velocity_y[j][i] = 0.0;
      g_velocity_x0[j][i] = 0.0;
      g_velocity_y0[j][i] = 0.0;
    }
  }
  // not finished yet!!!!!!!!!

  simulation(grid);
}

void Fluid_Grid::simulation(Fluid_Grid &grid) {
}