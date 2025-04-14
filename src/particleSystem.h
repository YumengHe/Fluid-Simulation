#ifndef PARTICLESYSTEM_H
#define PARTICLESYSTEM_H

#pragma once

#include <vector>
#include "particle.h"

class ParticleSystem {
public:
	// Simulation parameters (turnable for different fluid behaviors)
	float restDensity;			// p0: rest density for incompressibility
	float stiffness;			// k: pressure c
	float stiffnessNear;
	float viscosityLinear;
	float viscosityQuadratic;
	float springStiffness;
	float springYield;
	float plasticity;
	float radius;
	Vec2 gravity;

	std::vector<Particle> particles;
	std::vector<Spring> springs;

	ParticleSystem();
	void addParticle(const Vec2& pos, const Vec2& velocity = Vec2());
	void step(float dt);
	size_t getSpringCount() const;
};

#endif// PARTICLESYSTEM_H