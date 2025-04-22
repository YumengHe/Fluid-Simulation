#include "particle.h"

constexpr int DAM_PARTICLES = 10;
constexpr float VIEW_WIDTH = 800.f;
constexpr float VIEW_HEIGHT = 600.f;
constexpr float REST_DENS = 300.0f;     // rest density
constexpr float GAS_CONST = 2000.0f;    // const for equation of state
constexpr float H = 16.0f;              // smoothing length (kernel radius)
constexpr float VISC = 200.0f;          // viscosity constant
constexpr float DT = 0.0007f;           // integration timestep
constexpr float HSQ = H * H;            // radius^2 for optimization
constexpr float MASS = 2.5f;            // particle mass

// simulation parameters
constexpr float EPS = H;                // boundary epsilon
constexpr float BOUND_DAMPING = -0.5f;

const Vec2 G(0.0f, -9.81f); // External (gravitational) force

// Kernel constants for gradients
constexpr float POLY6 = 4.0f / (M_PI * (H * H * H * H * H * H * H * H));   // Poly6 kernel normalization
constexpr float SPIKY_GRAD = -10.0f / (M_PI * (H * H * H * H * H)); // Spiky gradient kernel
constexpr float VISC_LAP = 40.0f / (M_PI * (H * H * H * H * H));    // Viscosity Laplacian kernel

std::vector<Particle> particles;

// Dam break configuration for initialization of particles, can be switched out
void InitSPH() 
{
    std::cout << "Initializing dam break with " << DAM_PARTICLES << " particles" << std::endl;
    for (float y = EPS; y < VIEW_HEIGHT - EPS * 2.f; y += H) {
        for (float x = VIEW_WIDTH / 4; x <= VIEW_WIDTH / 2; x += H) {
            if (particles.size() < DAM_PARTICLES) {
                float jitter = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
                particles.push_back(Particle(x + jitter, y));
            }
            else{
                return;
            }
        }
    }
}

void ComputeDensityPressue()
{
    for (auto &pi : particles) {
        pi.rho = 0.f;
        for (auto &pj : particles) {    // summing density contibutions
            Vec2 rij = pj.x - pi.x;
            float r2 = rij.squaredNorm();

            if (r2 < HSQ) {
                pi.rho += MASS * POLY6 * pow(HSQ - r2, 3.f);
            }
        }
        pi.p = GAS_CONST * (pi.rho - REST_DENS);
    }
}

void ComputeForces()
{
    for (auto &pi : particles) {
        Vec2 fpress(0.f, 0.f);  // pressure force
        Vec2 fvisc(0.f, 0.f);   // viscosity force
        for (auto & pj : particles) {
            if (&pi == &pj) {
                continue;
            }

            Vec2 rij = pj.x - pi.x;
            float r = rij.norm();

            if (r < H) {
                // compute pressure force contribution
                fpress += -rij.normalized() * MASS * (pi.p + pj.p) / (2.f * pj.rho) * SPIKY_GRAD * pow(H - r, 3.f);
                // compute viscosity force contribution
                fvisc += VISC * MASS * (pj.v) / pj.rho * VISC_LAP * (H - r);
            }
        }
        Vec2 fgrav = G * MASS / pi.rho;
        pi.f = fpress + fvisc + fgrav;
    }
}

// Numerical Integration to update particle positions
void Integrate()
{
    for (auto &p : particles)
    {
        // forward Euler integration
        p.v += DT * p.f / p.rho;
        p.x += DT * p.v;

        // enforce boundary conditions
        if (p.x(0) - EPS < 0.f) {   // Left wall
            p.v(0) *= BOUND_DAMPING;
            p.x(0) = EPS;
        }
        if (p.x(0) + EPS > VIEW_WIDTH) {    // Right wall
            p.v(0) *= BOUND_DAMPING;
            p.x(0) = VIEW_WIDTH - EPS;
        } 
        if (p.x(0) - EPS < 0.f) {   // Bottom wall
            p.v(1) *= BOUND_DAMPING;
            p.x(1) = EPS;
        }
        if (p.x(1) + EPS > VIEW_HEIGHT) {    // Top wall
            p.v(1) *= BOUND_DAMPING;
            p.x(1) = VIEW_HEIGHT - EPS;
        } 
    }
}

void Update()
{
    ComputeDensityPressure();
    ComputeForces();
    Integrate();
}

// int main() {
//   // test for eigen
//   Eigen::MatrixXd m(2,2);
//   m(0,0) = 3;
//   m(1,0) = 2.5;
//   m(0,1) = -1;
//   m(1,1) = m(1,0) + m(0,1);
//   std::cout << m << std::endl;
// } 