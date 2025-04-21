#include <iostream>
#include <fstream>
#include <string>

class GridGenerator {
private:
    int g_width;
    int g_height;
    int g_num_iteration;
    float g_dt;
    float g_diffusion;
    float g_viscosity;
    float g_velocity_x;
    float g_velocity_y;
    float g_pressure;
    float g_density;

public:
    GridGenerator(int w = 5, int h = 5, int iter = 40, float dt = 0.1f,
                 float diff = 10.0f, float visc = 20.0f, float vel_x = 0.0f,
                 float vel_y = 0.0f, float pres = 1.0f, float dens = 2.0f)
        : g_width(w), g_height(h), g_num_iteration(iter), g_dt(dt),
          g_diffusion(diff), g_viscosity(visc), g_velocity_x(vel_x),
          g_velocity_y(vel_y), g_pressure(pres), g_density(dens) {}

    void generateGridFile(const std::string& filename) {
        std::ofstream file(filename);
        if (!file) {
            std::cerr << "Cannot create file: " << filename << std::endl;
            return;
        }

        // set precision and fixed notation
        file.precision(1);
        file << std::fixed;

        file << g_width << "\n";
        file << g_height << "\n";
        file << g_num_iteration << "\n";
        file << g_dt << "\n";
        file << g_diffusion << "\n";
        file << g_viscosity << "\n";

        for (int i = 0; i < g_width * g_height; ++i) {
            file << g_velocity_x << "\n";
        }

        for (int i = 0; i < g_width * g_height; ++i) {
            file << g_velocity_y << "\n";
        }

        for (int i = 0; i < g_width * g_height; ++i) {
            file << g_pressure << "\n";
        }

        for (int i = 0; i < g_width * g_height; ++i) {
            file << g_density << "\n";
        }

        file.close();
    }
};

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <filename>\n";
        return 1;
    }

    std::string filename = argv[1];
    
    // create a grid generator
    GridGenerator generator;
    generator.generateGridFile(filename);
    
    std::cout << "Created grid file: " << filename << std::endl;
    return 0;
}