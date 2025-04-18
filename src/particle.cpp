#include "particle.h"

Particle::Particle() : pos(), vel() {}
Particle::Particle(const Vec2& pos, const Vec2& vel) : pos(pos), vel(vel) {}

Spring::Spring(int i, int j, float restLen) : a(i), b(j), restLength(restLen) {}

// ------------------------------------------------------------
// Initialization & Simulation
// ------------------------------------------------------------
void Particle::initialization() {

}

void Particle::simulation() {
  
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