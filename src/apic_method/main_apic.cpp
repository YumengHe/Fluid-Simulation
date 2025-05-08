#include <GLUT/glut.h>
#include <chrono>
#include <iostream>

#include "array2_utils.h"
#include "fluidsim.h"

// Simulation Parameters
const int grid_resolution = 100;
const scalar domain_size = 100.0;
const scalar cfl_scale = 3.0;
const scalar frame_interval = 1.0 / 60.0;  // 60 FPS
const scalar max_step = 0.01;

FluidSim sim;
const Vector2s domain_center(50.0, 50.0);
const Vector2s origin(0.0, 0.0);

void setup_opengl() {
  glEnable(GL_DEPTH_TEST | GL_BLEND | GL_POINT_SMOOTH | GL_LINE_SMOOTH);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glHint(GL_POINT_SMOOTH_HINT | GL_LINE_SMOOTH_HINT, GL_NICEST);
  glClearColor(0.0, 0.0, 0.0, 1.0);
}

void render_scene() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0.0, domain_size, 0.0, domain_size);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  sim.render();
  glutSwapBuffers();
}

void handle_resize(int w, int h) { glViewport(0, 0, w, h); }

void key_callback(unsigned char key, int x, int y) {
  if (key == 27) std::exit(0);  // ESC
}

void tick(int) {
  glutPostRedisplay();
  glutTimerFunc(static_cast<int>(frame_interval * 1000.0), tick, 0);

  scalar dt = std::min(max_step, sim.compute_cfl() * cfl_scale);
  int steps = static_cast<int>(std::ceil(frame_interval / dt));
  dt = frame_interval / steps;

  for (int i = 0; i < steps; ++i) sim.advance(dt);
}

int main(int argc, char** argv) {
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(800, 800);
  glutCreateWindow("APIC Fluid Simulation");

  setup_opengl();
  sim.initialize(origin, domain_size, grid_resolution, grid_resolution, 1.0);
  sim.set_root_boundary(FluidSim::Boundary(domain_center, Vector2s(40.0, 40.0), FluidSim::BT_BOX, true));
  sim.update_boundary();
  sim.init_random_particles();

  glutDisplayFunc(render_scene);
  glutReshapeFunc(handle_resize);
  glutKeyboardFunc(key_callback);
  glutTimerFunc(1000, tick, 0);
  glutMainLoop();

  return 0;
}