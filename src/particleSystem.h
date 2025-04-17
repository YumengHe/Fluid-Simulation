#ifndef PARTICLESYSTEM_H
#define PARTICLESYSTEM_H

#pragma once

#include <vector>
#include "particle.h"

class ParticleSystem {
public:
	// Simulation parameters (turnable for different fluid behaviors)
	float restDensity;			// p0: rest density for incompressibility
	float stiffness;			// k: pressure constant for density
	float stiffnessNear;		// kNear: pressure constant for near density
	float viscosityLinear;		// σ: linear viscosity coefficien
	float viscosityQuadratic;	// β: quadratic viscosity coefficient
	float springStiffness;		// spring elasticity stiffness factor
	float springYield;			// μ: yield ratio (fraction of deformation tolerated without plastic change)
	float plasticity;			// α: plasticity constant (rate of rest-length change per timestep)
	float radius;				// h: interaction radius for neighbors
	Vec2 gravity;				// gravity acceleration vector

	std::vector<Particle> particles;
	std::vector<Spring> springs;

	ParticleSystem();
	void addParticle(const Vec2& pos, const Vec2& velocity = Vec2());
	void step(float dt);
	size_t getSpringCount() const;
};

#endif// PARTICLESYSTEM_H