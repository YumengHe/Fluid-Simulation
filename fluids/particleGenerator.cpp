#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include "../src/constants.h"

/**
 * @brief Generate particles in a grid pattern
 * @param rows Number of rows in the particle grid
 * @param cols Number of columns in the particle grid
 * @param spacing Distance between adjacent particles
 * @param outputFile Output file path for the particle data
 */
void generateParticles(
    int rows,
    int cols,
    float spacing,
    const std::string& outputFile = "particles.par"
) {
    // Open output file stream
    std::ofstream file(outputFile);
    
    // Write total number of particles as the first line
    int totalParticles = rows * cols;
    // file << totalParticles << "\n";
    
    // Calculate starting position (centered in view)
    float startX = (VIEW_WIDTH - (cols - 1) * spacing) / 2.0f;
    float startY = (VIEW_HEIGHT - (rows - 1) * spacing) / 2.0f;
    
    // Generate particles in a grid pattern
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            // Calculate particle position
            float x = startX + j * spacing;
            float y = startY + i * spacing;
            
            // Write position to file
            file << x << " " << y;
            // Only add newline if this is not the last particle
            if (!(i == rows-1 && j == cols-1)) {
                file << "\n";
            }
        }
    }
    
    file.close();
    
    // Print generation summary
    std::cout << "✅ Generated " << totalParticles << " particles to " << outputFile << std::endl;
}

int main() {
    // Example usage: Generate a 5x5 particle grid with 40 units spacing
    generateParticles(
        20,          // Number of rows
        20,          // Number of columns
        20.0,        // Spacing between particles
        "400particles.par"  // Output file name
    );
    return 0;
}