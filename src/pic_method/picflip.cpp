#include <GLUT/glut.h>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include "picflip.h"     
#include "constants.h"   // defines VIEW_WIDTH / VIEW_HEIGHT

int grid_size = 64;
float dt = 0.04f;
int pressure_iter = 40;
float flip_ratio = 0.96f;
bool simulation_start = false;

struct Vec2 {
    float x;
    float y;

    // Constructor
    Vec2(float x_val = 0.0f, float y_val = 0.0f) {
        x = x_val;
        y = y_val;
    }
    // Add
    Vec2 operator+(const Vec2& other) const {
        return Vec2(this->x + other.x, this->y + other.y);
    }
    // Subtract
    Vec2 operator-(const Vec2& other) const {
        return Vec2(this->x - other.x, this->y - other.y);
    }
    // Multiply 
    Vec2 operator*(float scalar) const {
        return Vec2(this->x * scalar, this->y * scalar);
    }
    Vec2& operator+=(const Vec2& other) {
        this->x += other.x;
        this->y += other.y;
        return *this;
    }
};


struct Particle {
    Vec2 pos; //Current position
    Vec2 vel; //Current velocity
};

struct Grid {
    // current grid velocity
    std::vector<std::vector<Vec2>> velocity;
    // old velocity
    std::vector<std::vector<Vec2>> velocity_old;
    // per cell mass
    std::vector<std::vector<float>> mass;
    std::vector<std::vector<float>> pressure;

    // Constructor
    Grid(int size) {
        velocity = std::vector<std::vector<Vec2>>(size, std::vector<Vec2>(size, Vec2(0, 0)));
        velocity_old = velocity;
        mass = std::vector<std::vector<float>>(size, std::vector<float>(size, 0.0f));
        pressure = mass;
    }

    // Reset grid
    void clear() {
        for (int i = 0; i < grid_size; ++i) {
            for (int j = 0; j < grid_size; ++j) {
                velocity[i][j] = Vec2(0, 0);
                mass[i][j]     = 0.0f;
                pressure[i][j] = 0.0f;
            }
        }
    }
};

std::vector<Particle> particles;
Grid grid(grid_size);

// Check if the grid cell [i, j] is fluid
bool isFluid(int i, int j) {
    if (i < 0 || i >= grid_size){
        return false;
    }
    if (j < 0 || j >= grid_size){
        return false;
    }
    return grid.mass[i][j] > 0.0f;
}

// Computes the B-spline weight
float bspline_weight(float x) {
    x = std::fabs(x);  // ensure x is non-negative
    if (x < 0.5f) {
        return 0.75f - x * x;
    } 
    else if (x < 1.5f) {
        float diff = 1.5f - x;
        return 0.5f * diff * diff;
    } 
    else {
        return 0.0f;
    }
}

// Initialize particles in a square block in the bottom-left of the grid
void initialize_particles() {
    particles.clear();  // clear any existing particles

    const float spacing = 0.3f; // Distance between particles
    float particle_area = grid_size * 0.5f;  // set particle area

    for (float x = 1.0f; x < particle_area; x += spacing) {
        for (float y = 1.0f; y < particle_area; y += spacing) {
            Vec2 position(x, y);
            Vec2 velocity(0.0f, 0.0f);
            particles.push_back(Particle{position, velocity});
        }
    }
}

// Transfer particle velocities to the grid using B-spline weights
void particle_to_grid() {
    grid.clear();  // Reset grid
    for (auto& p : particles) {
        // Find the base index
        int base_x = static_cast<int>(p.pos.x - 0.5f);
        int base_y = static_cast<int>(p.pos.y - 0.5f);
        // Loop 3x3 neighbor
        for (int i = -1; i <= 1; ++i) {
            for (int j = -1; j <= 1; ++j) {
                int grid_x = base_x + i;
                int grid_y = base_y + j;
                if (grid_x >= 0 && grid_x < grid_size && grid_y >= 0 && grid_y < grid_size) {
                    //B-spline weights for interpolation in x and y directions
                    float w_x = bspline_weight(p.pos.x - grid_x);
                    float w_y = bspline_weight(p.pos.y - grid_y);
                    float weight = w_x * w_y; //// Final weight
                
                    grid.velocity[grid_x][grid_y] += p.vel * weight;
                    grid.mass[grid_x][grid_y] += weight;
                }
            }
        }
    }

    // Normalize grid velocity by mass
    for (int i = 0; i < grid_size; ++i) {
        for (int j = 0; j < grid_size; ++j) {
            if (grid.mass[i][j] > 0.0f) {
                grid.velocity[i][j] = grid.velocity[i][j] * (1.0f / grid.mass[i][j]);
            }
        }
    }
}

// Apply forces (gravity) to grid velocities
void apply_grid_forces() {
    for (int i = 0; i < grid_size; ++i)
        for (int j = 0; j < grid_size; ++j)
            if (grid.mass[i][j] > 0)
                grid.velocity[i][j].y -= 9.8f * dt;
}


// Compute the velocity divergence at grid cell[i, j]
float compute_divergence(int i, int j) {
    float divergence = 0.0f;
    if (i + 1 < grid_size && isFluid(i + 1, j)) {
        divergence += grid.velocity[i + 1][j].x;  // right
    }
    if (i - 1 >= 0 && isFluid(i - 1, j)) {
        divergence -= grid.velocity[i - 1][j].x;  // left
    }
    if (j + 1 < grid_size && isFluid(i, j + 1)) {
        divergence += grid.velocity[i][j + 1].y;  // top
    }
    if (j - 1 >= 0 && isFluid(i, j - 1)) {
        divergence -= grid.velocity[i][j - 1].y;  // bottom
    }
    // central difference approximation
    return divergence * 0.5f;
}

// Solve for pressure and subtract its gradient to make the fluid incompressible.
void project_pressure() {
    // Solve ∇²p = divergence using Jacobi iterations
    for (int k = 0; k < pressure_iter; ++k) {
        std::vector<std::vector<float>> pressure_new = grid.pressure;
        for (int i = 1; i < grid_size - 1; ++i) {
            for (int j = 1; j < grid_size - 1; ++j) {
                if (isFluid(i, j)) {
                    float div = compute_divergence(i, j);  // Current cell divergence
                    int fluid_neighbors = 0;
                    float sum_pressure = 0.0f;

                    // Sum up pressure
                    if (isFluid(i + 1, j)) {
                        sum_pressure += grid.pressure[i + 1][j];
                        fluid_neighbors++;
                    }
                    if (isFluid(i - 1, j)) {
                        sum_pressure += grid.pressure[i - 1][j];
                        fluid_neighbors++;
                    }
                    if (isFluid(i, j + 1)) {
                        sum_pressure += grid.pressure[i][j + 1];
                        fluid_neighbors++;
                    }
                    if (isFluid(i, j - 1)) {
                        sum_pressure += grid.pressure[i][j - 1];
                        fluid_neighbors++;
                    }

                    // Update pressure
                    if (fluid_neighbors > 0) {
                        pressure_new[i][j] = (sum_pressure - div) / fluid_neighbors;
                    }
                }
            }
        }
        grid.pressure = pressure_new;
    }
    // Subtract pressure gradient from velocity
    for (int i = 1; i < grid_size - 1; ++i) {
        for (int j = 1; j < grid_size - 1; ++j) {
            if (isFluid(i, j)) {
                grid.velocity[i][j].x -= 0.5f * (grid.pressure[i + 1][j] - grid.pressure[i - 1][j]);
                grid.velocity[i][j].y -= 0.5f * (grid.pressure[i][j + 1] - grid.pressure[i][j - 1]);
            }
        }
    }
    // zero velocity at edges
    for (int i = 0; i < grid_size; ++i) {
        grid.velocity[i][0] = Vec2(0.0f, 0.0f);
        grid.velocity[i][grid_size - 1] = Vec2(0.0f, 0.0f);
        grid.velocity[0][i] = Vec2(0.0f, 0.0f);
        grid.velocity[grid_size - 1][i] = Vec2(0.0f, 0.0f);
    }
}

// Transfer velocity from the grid to particle
// Transfer velocity from the grid back to each particle using B-spline weights.
// Combines PIC and FLIP updates using a blending ratio.
void grid_to_particle() {
    for (Particle& p : particles) {
        int base_x = static_cast<int>(p.pos.x - 0.5f);
        int base_y = static_cast<int>(p.pos.y - 0.5f);
        Vec2 new_velocity(0.0f, 0.0f);
        Vec2 old_velocity(0.0f, 0.0f);
        // Loop over 3x3 neighboring cells
        for (int i = -1; i <= 1; ++i) {
            for (int j = -1; j <= 1; ++j) {
                int grid_x = base_x + i;
                int grid_y = base_y + j;
                if (grid_x >= 0 && grid_x < grid_size && grid_y >= 0 && grid_y < grid_size) {
                    float weight_x = bspline_weight(p.pos.x - grid_x);
                    float weight_y = bspline_weight(p.pos.y - grid_y);
                    float weight = weight_x * weight_y;

                    // Accumulate weighted velocity from grid
                    new_velocity += grid.velocity[grid_x][grid_y] * weight;
                    old_velocity += grid.velocity_old[grid_x][grid_y] * weight;
                }
            }
        }
        Vec2 flip_delta = new_velocity - old_velocity; //change in grid velocity
        // Blend FLIP and PIC
        p.vel = p.vel + flip_delta * flip_ratio + new_velocity * (1.0f - flip_ratio);
    }
}

// Move each particle forward
void advect_particles() {
    for (Particle& p : particles) {
        p.pos += p.vel * dt;// Update particle position
        // keep particles inside the grid
        if (p.pos.x < 1) {
            p.pos.x = 1;
            p.vel.x = 0;
        }
        if (p.pos.x > grid_size - 2) {
            p.pos.x = grid_size - 2;
            p.vel.x = 0;
        }
        if (p.pos.y < 1) {
            p.pos.y = 1;
            p.vel.y = 0;
        }
        if (p.pos.y > grid_size - 2) {
            p.pos.y = grid_size - 2;
            p.vel.y = 0;
        }
    }
}

// Perform simulation step
void step() {
    particle_to_grid();
    //Save current grid velocity for FLIP
    grid.velocity_old = grid.velocity;
    apply_grid_forces();
    project_pressure();
    grid_to_particle();
    advect_particles();
}

// Render all particles
void display() {
    glClear(GL_COLOR_BUFFER_BIT);
    glPointSize(2.0f);                    
    glBegin(GL_POINTS);
    glColor3f(0.2f, 0.6f, 1.0f); // particle color:blue
    for (Particle& p : particles) {
        // Normalize particle position to [-1, 1] OpenGL screen space
        float px = p.pos.x / grid_size * 2.0f - 1.0f;
        float py = p.pos.y / grid_size * 2.0f - 1.0f;
        glVertex2f(px, py);
    }

    glEnd();
    glutSwapBuffers();
}

void idle() {
    if (simulation_start) {
        step();
        glutPostRedisplay();
    }
}

// Handle keyboard input
void keyboard(unsigned char key, int x, int y) {
    if (key == 27) {
        exit(0);
    }
    if (key == 'p' || key == 'P') {
        simulation_start = true;
    }
}

void run_picflip_demo() {
    int argc = 1;
    char* argv[] = {(char*)"run"};

    initialize_particles();
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize((int)VIEW_WIDTH, (int)VIEW_HEIGHT);
    glutCreateWindow("PIC/FLIP Demo");
    glClearColor(0, 0, 0, 1);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glEnable(GL_POINT_SMOOTH);
    glutDisplayFunc(display);
    glutIdleFunc(idle);
    glutKeyboardFunc(keyboard);
    glutMainLoop();
}
