#include "particle.h"
#include "constants.h"

// solver parameters
const Vector2d G(0.0f, -9.81f); // External (gravitational) force
const static float REST_DENS = 300.0f; // rest density
constexpr float GAS_CONST = 2000.0f; // const for equation of state
constexpr float H = 16.0f;  // smoothing length (kernel radius)
constexpr float HSQ = H * H; // radius^2 for optimization
constexpr float MASS = 2.5f;            // particle mass
constexpr float VISC = 200.0f;          // viscosity constant
constexpr float DT = 0.0007f;           // integration timestep

// Kernel constants for gradients
constexpr float POLY6 = 4.0f / (M_PI * (H * H * H * H * H * H * H * H));   // Poly6 kernel normalization
constexpr float SPIKY_GRAD = -10.0f / (M_PI * (H * H * H * H * H)); // Spiky gradient kernel
constexpr float VISC_LAP = 40.0f / (M_PI * (H * H * H * H * H));    // Viscosity Laplacian kernel

// simulation parameters
constexpr float EPS = H;                // boundary epsilon
constexpr float BOUND_DAMPING = -0.5f;

// solver data
std::vector<Particle> particles;

// interaction
constexpr int DAM_PARTICLES = 500;

// --------------------------------------------------
// Initialization
// --------------------------------------------------
// Dam break configuration for initialization of particles, can be switched out
void InitSPH() {
    std::cout << "Initializing dam break with " << DAM_PARTICLES << " particles" << std::endl;
    for (float y = EPS; y < VIEW_HEIGHT - EPS * 2.f; y += H) {
        for (float x = VIEW_WIDTH / 4; x <= VIEW_WIDTH / 2; x += H) {
            if (particles.size() < DAM_PARTICLES) {
                float jitter = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
                particles.push_back(Particle(x + jitter, y));
            }
            else {
                return;
            }
        }
    }
}

// --------------------------------------------------
// Density and Pressure
// --------------------------------------------------
void ComputeDensityPressure() {
    for (auto &pi : particles) {
        pi.rho = 0.f;
        for (auto &pj : particles) {    // summing density contibutions
            Vector2d rij = pj.x - pi.x;
            float r2 = rij.squaredNorm();

            if (r2 < HSQ) {
                pi.rho += MASS * POLY6 * pow(HSQ - r2, 3.f);
            }
        }
        pi.p = GAS_CONST * (pi.rho - REST_DENS);
    }
}

// --------------------------------------------------
// Calculating Forces
// --------------------------------------------------
void ComputeForces() {
    for (auto &pi : particles) {
        Vector2d fpress(0.f, 0.f);  // pressure force
        Vector2d fvisc(0.f, 0.f);   // viscosity force
        for (auto & pj : particles) {
            if (&pi == &pj) {
                continue;
            }

            Vector2d rij = pj.x - pi.x;
            float r = rij.norm();

            if (r < H) {
                // compute pressure force contribution
                fpress += -rij.normalized() * MASS * (pi.p + pj.p) / (2.f * pj.rho) * SPIKY_GRAD * pow(H - r, 3.f);
                // compute viscosity force contribution
                fvisc += VISC * MASS * (pj.v) / pj.rho * VISC_LAP * (H - r);
            }
        }
        Vector2d fgrav = G * MASS / pi.rho;
        pi.f = fpress + fvisc + fgrav;
    }
}

// --------------------------------------------------
// Numerical Integration
// --------------------------------------------------
// Numerical Integration to update particle positions
void Integrate() {
    for (auto &p : particles) {
        // forward Euler integration
        p.v += DT * p.f / p.rho;
        p.x += DT * p.v;

        // enforce boundary conditions
        if (p.x(0) - EPS < 0.f) { // Left wall
            p.v(0) *= BOUND_DAMPING;
            p.x(0) = EPS;
        }
        if (p.x(0) + EPS > VIEW_WIDTH) { // Right wall
            p.v(0) *= BOUND_DAMPING;
            p.x(0) = VIEW_WIDTH - EPS;
        } 
        if (p.x(1) - EPS < 0.f) { // Bottom wall
            p.v(1) *= BOUND_DAMPING;
            p.x(1) = EPS;
        }
        if (p.x(1) + EPS > VIEW_HEIGHT) { // Top wall
            p.v(1) *= BOUND_DAMPING;
            p.x(1) = VIEW_HEIGHT - EPS;
        } 
    }
}

// --------------------------------------------------
// Solver
// --------------------------------------------------
void Update() {
    ComputeDensityPressure();
    ComputeForces();
    Integrate();
    std::cout << "iterating..." << std::endl;
}
