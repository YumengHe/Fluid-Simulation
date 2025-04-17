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
void diffuse(int N, int b, std::vector< std::vector<float> > &x, std::vector< std::vector<float> > &x0, float diff, float dt, int num_iteration) {
  float a = dt * diff * N * N;
  lin_solve(b, x, x0, a, 1 + 4 * a, num_iteration, N);
}

// ------------------------------------------------------------
// Projection
// ------------------------------------------------------------
void project(int N, std::vector< std::vector<float> > &velocity_x, std::vector< std::vector<float> > &velocity_y, std::vector< std::vector<float> > &p, std::vector< std::vector<float> > &div, int num_iteration) {
  float h = 1.0 / N;
  for (int i = 1; i <= N; i ++) {
    for (int j = 1; j <= N; j ++) {
      div[i][j] = -0.5 * h * (velocity_x[i + 1][j] - velocity_x[i - 1][j] + velocity_y[i][j + 1] - velocity_y[i][j - 1]);
      p[i][j] = 0;
    }
  }
  set_bnd(N, 0, div);
  set_bnd(N, 0, p);

  lin_solve(0, p, div, 1, 4, num_iteration, N);

  for (int i = 1; i <= N; i ++) {
    for (int j = 1; j <= N; j ++) {
      velocity_x[i][j] -= 0.5 * (p[i + 1][j] - p[i - 1][j]) / h;
      velocity_y[i][j] -= 0.5 * (p[i][j + 1] - p[i][j - 1]) / h;
    }
  }
  set_bnd(N, 1, velocity_x);
  set_bnd(N, 2, velocity_y);
}

// ------------------------------------------------------------
// Helper functions
// ------------------------------------------------------------
// apply boundary conditions (e.g. walls)
// N: number of grid
// b: 0->scalar, 1->horizontal velocity, 2->vertical velocity
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
// c: ??
// num_iteration: number of iterations (more iteration -> more accurate)
// N: number of grid
void lin_solve(int b, std::vector< std::vector<float> > &x, std::vector< std::vector<float> > &x0, float a, float c, int num_iteration, int N) {
  for (int k = 0; k < num_iteration; k ++) {
    for (int j = 1; j <= N; j ++) { // y
      for (int i = 1; i <= N; i++) { // x
        x[j][i] = (x0[j][i] + a * (x[j - 1][j] + x[j + 1][i] + x[j][i - 1] + x[j][i + 1])) / c;
      }
    }
    set_bnd(N, b, x);
  }
};

// for user input (such as mouse drag), not necessary for our program but keep it here for now
// N: number of grid
// s: source of input density
// x: density
void add_source(int N, std::vector< std::vector<float> > &x, std::vector< std::vector<float> > &s, float dt) {
  for (int i = 0; i < N + 2; i ++) {
    for (int j = 0; j < N + 2; j ++) {
      x[i][j] += dt * s[i][j];
    }
  }
}

// ------------------------------------------------------------
// Core functions
// ------------------------------------------------------------
void dens_step() {

}

void vel_step() {
  
}

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