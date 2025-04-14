#pragma once
#include <vector>
#include <include/eigen-3.4.0/Eigen/Dense>

class Fluid_Particle{
// Represents a fluid particle with position and velocity
struct Particle {
    Vec2 pos;
    Vec2 vel;
    // Constructors
    Particle();
    Particle(const Vec2& pos, const Vec2& vel = Vec2());
};

// Represents an elastic spring connecting two particles (indices a and b)
struct Spring {
    int a;          // index of first particle
	int b;		    // index of second particle
    float restLength;
	// Constructor
	Spring(int i, int j, float restLength);
};

};
