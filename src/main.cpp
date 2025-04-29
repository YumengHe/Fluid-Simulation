#include <stdio.h>     /* printf */
#include <GLUT/glut.h> /* glut graphics library */
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include "input.h"
#include "grid.h"
#include "particle.h"
#include "constants.h"
std::vector<Particle> particles;  // 声明这是一个外部变量
#include "./pic_method/apic.h"
#include "./pic_method/pic.h"
int mode = 0; // 0 grid, 1 particle, 2 PIC, 3 APIC


void display(){ // change to particles later
    glClear(GL_COLOR_BUFFER_BIT);
    if (mode == 0) {
        // Draw grid
        glBegin(GL_QUADS);
        float cell_width = 1.0f / current_grid->g_width;
        float cell_height = 1.0f / current_grid->g_height;
        // Commented out: Print grid size
        for (int i = 0; i < current_grid->g_height; i++) {
            for (int j = 0; j < current_grid->g_width; j++) {
                // Commented out: Print density values
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
            // Commented out: Line break
        }
        glEnd();
    } else if (mode == 1) {  // Check if there are particles to display
        // Enable point smoothing for circular particles
        glEnable(GL_POINT_SMOOTH);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
        glPointSize(8.0f);  // Set point size
        glBegin(GL_POINTS); 
        glColor3f(0.0f, 0.4f, 1.0f);  // Blue particles
        
        for(const auto& p : particles) { 
            // Map particle coordinates to [0,1] range
            float x = p.x[0] / VIEW_WIDTH; 
            float y = p.x[1] / VIEW_HEIGHT; 
            glVertex2f(x, y); 
        } 
        glEnd();
        
        // Disable point smoothing after drawing
        glDisable(GL_BLEND);
        glDisable(GL_POINT_SMOOTH);
    }else if (mode == 3) {  // 添加APIC粒子的渲染
        // Enable point smoothing for circular particles
        glEnable(GL_POINT_SMOOTH);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
        glPointSize(8.0f);
        glBegin(GL_POINTS); 
        glColor3f(0.0f, 0.6f, 0.8f);  // 使用不同的蓝色来区分APIC粒子
        
        for(const auto& p : apic_particles) { 
            float x = p.x; 
            float y = p.y; 
            glVertex2f(x, y); 
        } 
        glEnd();
        
        glDisable(GL_BLEND);
        glDisable(GL_POINT_SMOOTH);
    }
    else if (mode == 2) {  // 添加PIC粒子的渲染
        glEnable(GL_POINT_SMOOTH);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
        glPointSize(8.0f);
        glBegin(GL_POINTS); 
        glColor3f(1.0f, 0.4f, 0.0f);  // 使用橙色来区分PIC粒子
        
        for(const auto& p : pic_particles) { 
            float x = p.x; 
            float y = p.y; 
            glVertex2f(x, y); 
        } 
        glEnd();
        
        glDisable(GL_BLEND);
        glDisable(GL_POINT_SMOOTH);
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
    if (argc > 3){
        std::cerr << "Please input single particle/grid file" << std::endl;
        return 1;
    }

    if (argc == 1) {
        mode = 1;
        InitSPH();
    }
    else {
        std::string filename = argv[1];
        if (endsWith(filename, ".grid")){
            mode = 0;
            std::cout << "Loading a .grid file\n";
            static Fluid_Grid grid;
            loadGrid(filename, grid); 
            current_grid = &grid; 
        }
        else if (endsWith(filename, ".par")){
            mode = 1;
            std::cout << "Loading a .par file\n";
            loadParticles(filename);
        }
        else if (endsWith(filename, ".apic")){
            mode = 3;
            std::cout << "Loading a .apic file\n";
            loadAPIC(filename); // Load APIC data
        }
        else if (endsWith(filename, ".pic")) {
            mode = 2;
            std::cout << "Loading a .pic file\n";
            loadPIC(filename); // 加载PIC数据
        }
        else{
            std::cerr << "Unsupported file type. Please use a .grid or .par file.\n";
            return 1;
        }
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
