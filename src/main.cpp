#include <stdio.h>     /* printf */
#include <GLUT/glut.h> /* glut graphics library */
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include "input.h"
#include "grid.h"

Fluid_Grid* current_grid = nullptr; 

void display(){ // change to particles later
    glClear(GL_COLOR_BUFFER_BIT);
    if (current_grid) {
        // Draw grid
        glBegin(GL_QUADS);
        float cell_width = 1.0f / current_grid->g_width;
        float cell_height = 1.0f / current_grid->g_height;
        // std::cout << "Grid size: " << current_grid->g_width << "x" << current_grid->g_height << std::endl;  // 注释掉
        for (int i = 0; i < current_grid->g_height; i++) {
            for (int j = 0; j < current_grid->g_width; j++) {
                // std::cout << current_grid->g_density[i][j] << " ";  // 注释掉

                float x = j * cell_width;
                float y = (current_grid->g_height - 1 - i) * cell_height;
                
                float density = current_grid->g_density[i][j];
                float color = density / 2.0f;
                color = std::min(1.0f, std::max(0.0f, color));
                
                glColor3f(color, color, color);
                
                glVertex2f(x, y);
                glVertex2f(x + cell_width, y);
                glVertex2f(x + cell_width, y + cell_height);
                glVertex2f(x, y + cell_height);
            }
            // std::cout << std::endl;  // 注释掉
        }
        glEnd();
    }
    glutSwapBuffers();
}


void initGL(){
    glClearColor(0, 0, 0, 1); // black background
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, 1, 0, 1); // 2D coordinate system from (0,0) to (1,1)
}

int main(int argc, char **argv){
    if (argc < 2){
        std::cerr << "Please input initial particle/grid file" << std::endl;
        return 1;
    }

    std::string filename = argv[1];
    if (endsWith(filename, ".grid")){
        std::cout << "Loading a .grid file\n";
        static Fluid_Grid grid;
        loadGrid(filename, grid); 
        current_grid = &grid; 
    }
    else if (endsWith(filename, ".par")){
        std::cout << "Loading a .par file\n";
        // loadParticles(filename);
    }
    else{
        std::cerr << "Unsupported file type. Please use a .grid or .par file.\n";
        return 1;
    }

    glutInit(&argc, argv); // initialize GLUT
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB); // double buffering and RGB color
    glutInitWindowSize(800, 800); // window size
    glutCreateWindow("2D OpenGL Fluid Framework");

    initGL();
    glutDisplayFunc(display);
    glutIdleFunc(idle);
    glutMouseFunc(mouseClick);
    glutKeyboardFunc(handleKeypress);

    glutMainLoop();
    return 0;
}
