#pragma once

#include <vector>
#include <iostream>
#include "../include/Eigen/Dense"

using namespace Eigen;

// Represents a fluid particle with position and velocity
struct Particle {
    Vector2d x; // position
    Vector2d v; // velocity
    Vector2d f; // force
    float rho; // density
    float p; // pressure

    // constructor
    Particle(float _x, float _y){
        x(_x, _y);
        v(0.0, 0.0);
        rho = 0.0; 
        p = 0.0;
    }
};

// physics functions
void InitSPH();
void Integrate();
void ComputeDensityPressure();
void ComputeForces();
void Update();
