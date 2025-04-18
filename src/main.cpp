#include <stdio.h>     /* printf */
#include <GLUT/glut.h> /* glut graphics library */
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include "input.h"
#include "../include/json.hpp"

using json = nlohmann::json;

struct Particle{
    float x, y;
};

std::vector<Particle> particles;

// read particles.json 
void loadParticles(const std::string &filename){
    std::ifstream file(filename);
    if (!file){
        std::cerr << "Failed to open " << filename << std::endl;
        return;
    }

    json data;
    file >> data;

    for (const auto &p : data["particles"]){
        Particle particle;
        particle.x = p["position"][0];
        particle.y = p["position"][1];
        particles.push_back(particle);
    }

    std::cout << "Loaded " << particles.size() << " particles.\n";
}

void display(){ // change to particles later
    glClear(GL_COLOR_BUFFER_BIT);

    glColor3f(0.0f, 0.6f, 1.0f); // blue color
    glPointSize(8.0f);
    glBegin(GL_POINTS);
    for (const auto &p : particles){
        glVertex2f(p.x, p.y); 
    }
    glEnd();

    glutSwapBuffers();
}

void initGL(){
    glClearColor(0, 0, 0, 1); // black background
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, 1, 0, 1); // 2D coordinate system from (0,0) to (1,1)
}

int main(int argc, char **argv){
    glutInit(&argc, argv); // initialize GLUT
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB); // double buffering and RGB color
    glutInitWindowSize(800, 800); // window size
    glutCreateWindow("2D OpenGL Fluid Framework");

    initGL();
    loadParticles("./particles.json");
    glutDisplayFunc(display);
    glutIdleFunc(idle);
    glutMouseFunc(mouseClick);
    glutKeyboardFunc(handleKeypress);

    glutMainLoop();
    return 0;
}
