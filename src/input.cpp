#include <GLUT/glut.h>
#include <cstdio>
#include <cstdlib>
#include "input.h"

void mouseClick(int button, int state, int x, int y)
{
    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
    {
        // convert screen coord (x, y) to OpenGL 0~1 region
        float xf = (float)x / glutGet(GLUT_WINDOW_WIDTH);
        float yf = 1.0f - (float)y / glutGet(GLUT_WINDOW_HEIGHT); // flip y coord

        printf("🖱️ Click at: screen = (%d, %d), normalized = (%.3f, %.3f)\n", x, y, xf, yf);
    }
}

void handleKeypress(unsigned char key, int x, int y)
{
    if (key == 27){            // 27 is the ASCII value of ESC
        exit(0); // esc
    }
}

void idle()
{
    glutPostRedisplay(); // request a redraw when cpu is idle
}

bool endsWith(const std::string &str, const std::string &suffix)
{
    return str.size() >= suffix.size() &&
           str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

// void loadGrid(const std::string &filename, Fluid_Grid &grid){

// }