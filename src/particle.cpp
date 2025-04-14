#include "Particle.h"

Particle::Particle() : pos(), vel() {}
Particle::Particle(const Vec2& pos, const Vec2& vel) : pos(pos), vel(vel) {}

Spring::Spring(int i, int j, float restLen) : a(i), b(j), restLength(restLen) {}