#include "particleSystem.h"
#include <cmath>

ParticleSystem::ParticleSystem() {
	restDensity = 5.0f;
	stiffness = 0.5f;
	stiffnessNear = 0.5f;
	viscosityLinear = 0.1f;
	viscosityQuadratic = 0.1f;
	springStiffness = 0.5f;
	springYield = 0.1f;
	plasticity = 0.3f;
	radius = 0.15f;
	gravity = Vec2(0.0f, -9.8f);
}

void ParticleSystem::addParticle(const Vec2& pos, const Vec2& vel) {
	particles.emplace_back(pos, vel);
}

size_t ParticleSystem::getSpringCount() const {
	return springs.size();
}

void ParticleSystem::step(float dt) {

}