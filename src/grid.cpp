#include "grid.h"
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
// Projection
// ------------------------------------------------------------
void Fluid_Grid::project() {
  
}

// ------------------------------------------------------------
// Helper functions
// ------------------------------------------------------------
// apply boundary conditions (e.g. walls)
// N: number of grid
// b: case
// x: velocity/density
void set_bnd(int N, int b, std::vector< std::vector<float> > &x) {
  for (int i = 0; i <= N; i ++) {
    if (b == 1) {
      x[0][i] = -x[1][i]; // upper wall
      x[N + 1][i] = -x[N][i]; // lower wall
    } else {
      x[0][i] = x[1][i]; // upper wall
      x[N + 1][i] = x[N][i]; // lower wall
    }

    if (b == 2) {
      x[i][0] = -x[i][1]; // left wall
      x[i][N + 1] = -x[i][N]; // right wall
    } else {
      x[i][0] = x[i][1]; // left wall
      x[i][N + 1] = x[i][N]; // right wall
    }

    // 4 corners
    x[0][0] = 0.5 * (x[1][0] + x[0][1]);
    x[0][N + 1] = 0.5 * (x[1][N + 1] + x[0][N]);
    x[N + 1][0] = 0.5 * (x[N][0] + x[N + 1][1]);
    x[N + 1][N + 1] = 0.5 * (x[N][N + 1] + x[N + 1][N]);
  }
}

// linear solver
// b: case
// x: velocity/density
// x0: old velocity/density
// a: diffusion rate
// num_iteration: number of iterations (more iteration -> more accurate)
// N: number of grid
void lin_solve(int b, std::vector< std::vector<float> > &x, std::vector< std::vector<float> > &x0, float a, int num_iteration, int N) {
  for (int k = 0; k < num_iteration; k ++) {
    for (int j = 1; j <= N; j ++) { // y
      for (int i = 1; i <= N; i++) { // x
        x[j][i] = (x0[j][i] + a * (x[j - 1][j] + x[j + 1][i] + x[j][i - 1] + x[j][i + 1])) / (1 + 4 * a);
      }
    }
    set_bnd(N, b, x);
  }
};

// ------------------------------------------------------------
// Initialization & Simulation
// ------------------------------------------------------------
void Fluid_Grid::initialization(){
  // note: we only consider equal length of each side -> width = height
  // +2 -> the grid contains an extra layer of cells to account for the boundary conditions
  g_width = 5 + 2;
  g_height = 5 + 2;

  g_dt = 0.1;
  g_diffusion = 10;
  g_viscosity = 10;
  g_num_iteration = 20;

  // resize
  // g_velocity[row/width][column/height]
  g_velocity_x.resize(g_height, std::vector<float>(g_width, 0.0f));
  g_velocity_y.resize(g_height, std::vector<float>(g_width, 0.0f));
  g_velocity_x0.resize(g_height, std::vector<float>(g_width, 0.0f));
  g_velocity_y0.resize(g_height, std::vector<float>(g_width, 0.0f));

  // initialize velocity
  float velocity = 0.0;
  for (int j = 0; j < g_height; j ++) {
    for (int i = 0; i < g_width; i ++) {
      g_velocity_x[j][i] = velocity;
      g_velocity_y[j][i] = velocity;
      g_velocity_x0[j][i] = velocity;
      g_velocity_y0[j][i] = velocity;
    }
  }
}

void Fluid_Grid::simulation() {
  std::cout << *this << std::endl;
}

// temporary main function
int main() {
  Fluid_Grid grid;
  grid.initialization();
  grid.simulation();
  return 0;
}