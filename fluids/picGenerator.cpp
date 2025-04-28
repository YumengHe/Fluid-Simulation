#include <iostream>
#include <fstream>
#include <vector>
#include "../src/constants.h"

/**
 * @brief Generate PIC particles in grid mode
 * @param dt Time step
 * @param grid_width Grid width
 * @param grid_height Grid height
 * @param rows Number of rows in the particle grid
 * @param cols Number of columns in the particle grid
 * @param spacing Distance between adjacent particles
 * @param outputFile Output file path for particle data
 */
void generatePICParticles(
    float dt,
    int grid_width,
    int grid_height,
    int rows,
    int cols,
    float spacing,
    const std::string& outputFile = "particles.pic"
) {
    std::ofstream file(outputFile);
    
    // Write time step
    file << dt << "\n";
    
    // Write grid information
    file << grid_width << " " << grid_height << "\n";
    
    // Calculate starting position (starting from bottom left corner)
    float grid_dx = VIEW_WIDTH / grid_width;
    float startX = grid_dx;  // Start from the first grid
    float startY = grid_dx;  // Start from the first grid
    
    // Generate particles in grid mode
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            // Calculate particle position
            float x = startX + j * spacing;
            float y = startY + (rows - 1 - i) * spacing;
            
            // Initial velocity set to 0
            float vx = 0.0f;
            float vy = 0.0f;
            
            // Write particle data: position, velocity
            file << x << " " << y << " "          // Position
                 << vx << " " << vy << "\n";      // Velocity
        }
    }
    
    file.close();
    
    // Print generation summary
    std::cout << "✅ PIC configuration file generated to " << outputFile << std::endl;
    std::cout << "   Grid size: " << grid_width << "x" << grid_height << std::endl;
    std::cout << "   Number of particles: " << rows * cols << std::endl;
}

int main() {
    // Example usage
    generatePICParticles(
        0.01f,      // dt: Time step
        50,         // number of grids in x direction
        50,         // number of grids in y direction
        5,          // rows: Number of particle rows
        5,          // cols: Number of particle columns
        20.0f,      // spacing: Particle spacing
        "test.pic"  // Output file name
    );
    return 0;
}