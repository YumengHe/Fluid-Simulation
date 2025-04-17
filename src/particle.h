#pragma once
#include <vector>
#include <iostream>
#include "Eigen/Dense"

using Vec2 = Eigen::Vector2f;

// Represents a fluid particle with position and velocity
struct Particle {
    // physics parameter
    Vec2 pos;
    Vec2 vel;
    Vec2 force;
    float density;
    float pressure;

    // constructors
    Particle();
    Particle(const Vec2& pos, const Vec2& vel = Vec2());

    // initialize
    void InitSPH();
};

// physics functions
void ComputeDensityPressure();
void ComputeForces();
void Integrate();

// core functions
void Update();

// Represents an elastic spring connecting two particles (indices a and b)
struct Spring {
    int a;          // index of first particle
	int b;		    // index of second particle
    float restLength;
	// Constructor
	Spring(int i, int j, float restLength);
};
