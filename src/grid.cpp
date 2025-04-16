#include "grid.h"
#include <cmath>

// ------------------------------------------------------------
// Constructor
// ------------------------------------------------------------
Fluid_Grid::Fluid_Grid() {

}

// looks like not necessary, but keep it here for now
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
// N: number of grid
// b: case
// d: destination
// d0: source
// v_x: velocity_x
// v_y: velocity_y
// dt: time step
void advect(int N, int b, std::vector< std::vector<float> > d, std::vector< std::vector<float> > d0, std::vector< std::vector<float> > v_x, std::vector< std::vector<float> > v_y, float dt) {
  float dt0 = dt * N;

  for (int i = 1; i <= N; i ++) {
    for (int j = 1; j <= N; j ++) {
      // backward trace step
      float x = i - dt0 * v_x[i][j];
      float y = j - dt0 * v_y[i][j];

      if (x < 0.5) x = 0.5; // clamp
      if (x > N + 0.5) x = N + 0.5; // clamp
      int i0 = (int)x;
      int i1 = i0 + 1;

      if (y < 0.5) y = 0.5; // clamp
      if (y > N + 0.5) y = N + 0.5; // clamp
      int j0 = (int)y;
      int j1 = j0 + 1;

      float s1 = x - i0;
      float s0 = 1 - s1;
      float t1 = y - j0;
      float t0 = 1 - t1;
      // bilinear interpolation
      d[i][j] = s0 * (t0 * d0[i0][j0] + t1 * d0[i0][j1]) + s1 * (t0 * d0[i1][j0] + t1 * d0[i1][j1]);
    }
  }
  set_bnd(N, b, d);
}

// ------------------------------------------------------------
// Diffusion
// ------------------------------------------------------------
// N: number of grid
// b: case
// x: destination
// x0: source
// diff: diffusion coefficient
// dt: time step
void diffuse(int N, int b, std::vector< std::vector<float> > &x, std::vector< std::vector<float> > &x0, float diff, float dt) {
  float a = dt * diff * N * N;
  set_bnd(N, b, x);
}

// ------------------------------------------------------------
// Projection
// ------------------------------------------------------------
void project() {
  
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
  int N = 5;
  g_width = N + 2;
  g_height = N + 2;

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