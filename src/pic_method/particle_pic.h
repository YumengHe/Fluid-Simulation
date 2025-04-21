#pragma once
struct Particle_PIC{
    float x; //The position in the x axis
    float y; //The position in the y axis

    float velocity_x;//The velocity in the x direction
    float velocity_y;//The velocity in the y direction
    
    //constructor
    Particle_PIC(float position_x,float position_y,float v_x,float v_y){
        x=position_x;
        y=position_y;
        velocity_x=v_x;
        velocity_y=v_y;
    }
};