#pragma once
#include <vector>
#include <iostream>
#include "../include/Eigen/Dense"

using Vec2 = Eigen::Vector2f;

// Represents a fluid particle with position and velocity
struct Particle {
    Particle(float _x, float _y) : x(_x, _y), v(0.f, 0.f), rho(0), p(0.f) {}
    Vec2 x, v, f;   // position, velocity, force
    float rho, p;   // density, pressure
};

// physics functions
void InitSPH();
void Integrate();
void ComputeDensityPressure();
void ComputeForces();
void Update();
