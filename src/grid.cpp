#include <grid.h>
#include <cmath>

// ------------------------------------------------------------
// Constructor
// ------------------------------------------------------------
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