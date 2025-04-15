
/*
Linux c console program
gcc f.c -lglut -lGL
./a.out
*/
// int main(int argc, char **argv)
// {
//     glutInit(&argc, argv);
//     glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
//     glutCreateWindow("red 3D lighted cube");
//     printf("GL_VERSION = %s\n", glGetString(GL_VERSION)); /* GL_VERSION = 2.1.2 NVIDIA 195.36.24 */
//     return 0;
// }
// main.cpp
#include <stdio.h>     /* printf */
#include <GLUT/glut.h> /* glut graphics library */
#include <cstdlib>

void display() // change to particles later
{
    glClear(GL_COLOR_BUFFER_BIT);

    // example: draw a test triangle
    glColor3f(0.2f, 0.7f, 1.0f);
    glBegin(GL_TRIANGLES);
    glVertex2f(0.25f, 0.25f);
    glVertex2f(0.75f, 0.25f);
    glVertex2f(0.5f, 0.75f);
    glEnd();

    glutSwapBuffers();
}

void idle()
{
    glutPostRedisplay(); // request a redraw when cpu is idle
}

void initGL()
{
    glClearColor(0, 0, 0, 1); // black background
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, 1, 0, 1); // 2D coordinate system from (0,0) to (1,1)
}

int main(int argc, char **argv)
{
    glutInit(&argc, argv); // initialize GLUT
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB); // double buffering and RGB color
    glutInitWindowSize(800, 800); // window size
    glutCreateWindow("2D OpenGL Fluid Framework");

    initGL();

    glutDisplayFunc(display);
    glutIdleFunc(idle);

    glutMainLoop();
    return 0;
}
