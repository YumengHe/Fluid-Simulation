#pragma once
#include <iostream>
#include "grid.h"

extern Fluid_Grid* current_grid;

void mouseClick(int button, int state, int x, int y);
void handleKeypress(unsigned char key, int x, int y);
void idle();
bool endsWith(const std::string &str, const std::string &suffix);
void loadParticles(const std::string &filename);
void loadGrid(const std::string &filename, Fluid_Grid &grid);
void loadPIC(const std::string &filename); 