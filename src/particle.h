#pragma once

#include <vector>
#include <iostream>
#include <GLUT/glut.h>
#include "../include/Eigen/Dense"

using namespace Eigen;
using namespace std;

// Represents a fluid particle with position and velocity
struct Particle {
    Vector2d x; // position
    Vector2d v; // velocity
    Vector2d f; // force
    float rho; // density
    float p; // pressure

    // constructor
    Particle(float _x, float _y) : x(_x, _y), v(0.f, 0.f), f(0.f, 0.f), rho(0), p(0.f) {}
};

// physics functions
void InitSPH(void);
void Integrate(void);
void ComputeDensityPressure(void);
void ComputeForces(void);
void Update(void);

extern std::vector<Particle> particles;
