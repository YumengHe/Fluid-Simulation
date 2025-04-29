#pragma once
#include "../../include/Eigen/Dense"

struct Particle_PIC {
    float x; //The position in the x axis
    float y; //The position in the y axis

    float velocity_x;//The velocity in the x direction
    float velocity_y;//The velocity in the y direction
    float new_velocity_x, new_velocity_y;
    // APIC affine velocity field
    Eigen::Matrix2f B; // affine matrix for velocity field
    
    //constructor
    Particle_PIC(float x = 0.0f, float y = 0.0f, float vx = 0.0f, float vy = 0.0f)
        : x(x), y(y), velocity_x(vx), velocity_y(vy), new_velocity_x(0.0f), new_velocity_y(0.0f) {
        B.setZero();
    }
};