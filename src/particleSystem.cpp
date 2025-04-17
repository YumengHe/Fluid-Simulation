#include "particleSystem.h"
#include "particle.h"
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
	size_t n = particles.size();
	if (n == 0) return;

	// Apply gravity to all particles
	// Algorithm 1 (3)
	for (size_t i = 0; i < n; i++) {
		particles[i].vel += gravity * dt;
	}

	// Apply pairwise viscosity damping
	// Algorithm 5 (5.3)
	for (size_t i = 0; i < n; i++) {
		for (size_t j = i + 1; j < n; j++) {
			Vec2 rij = particles[j].pos - particles[i].pos;	// relative position
			float dist2 = rij.squaredNorm();	// squared distance
			if (dist2 < radius * radius) {		// if within interaction radius
				float dist = std::sqrt(dist2);
				if (dist < 1e-6f) continue;
				Vec2 dir = rij / dist;	// direction vector
				float relativeVel = (particles[j].vel - particles[i].vel).dot(dir);	// relative velocity along the direction
				if (relativeVel > 0) {
					float damping = (1 - dist / radius) * (viscosityLinear * relativeVel + viscosityQuadratic * relativeVel * relativeVel);
					Vec2 impulse = dir * (damping * dt);
					particles[i].vel += impulse * 0.5f;
					particles[j].vel -= impulse * 0.5f;
				}
			}
		}
	}

	// Euler prediction step
	// Algorithm 1 (3)
	std::vector<Vec2> prevPos(n);
	for (size_t i = 0; i < n; i++) {
		prevPos[i] = particles[i].pos;				// save previous position
		particles[i].pos += particles[i].vel * dt;	// advance to predicted position
	}

	// Viscoelastic plasticity control, spring adjustment
	// Algorithm 4 (5.2)
	for (size_t i = 0; i < n; i++) {
		for (size_t j = i + 1; j < n; j++) {	// for each neighbor pair ij, (i < j)
			Vec2 rij = particles[j].pos - particles[i].pos;
			if (rij.squaredNorm() < radius * radius) {
				bool exists = false;
				for (const Spring& s : springs) {
					if ((s.a == (int)i && s.b == (int)j) || (s.a == (int)j && s.b == (int)i)) {
						exists = true;
						break;
					}
				}
				if (!exists)
					springs.emplace_back((int)i, (int)j, radius);
			}
		}
	}

	for (auto it = springs.begin(); it != springs.end();) {
		int i = it->a, j = it->b;
		Vec2 rij = particles[j].pos - particles[i].pos;
		float dist = rij.norm();
		float L0 = it->restLength;
		float yieldThresh = springYield * L0;
		if (dist > L0 + yieldThresh)
			it->restLength += plasticity * (dist - L0 - yieldThresh);
		else if (dist < L0 - yieldThresh)
			it->restLength -= plasticity * (L0 - yieldThresh - dist);

		if (it->restLength > radius) {
			it = springs.erase(it);
			continue;
		}
		++it;
	}

	// Spring elasticity: particle displacement from spring stretch
	// Algorithm 3 (5.1)
	for (const Spring& s : springs) {	// for each spring
		int i = s.a, j = s.b;
		Vec2 rij = particles[j].pos - particles[i].pos;
		float dist = rij.norm();
		if (dist < 1e-6f) continue;
		float stretch = 1.0f - (s.restLength / dist);
		Vec2 springDisp = rij * (springStiffness * stretch * dt * dt);
		particles[i].pos += springDisp * 0.5f;
		particles[j].pos -= springDisp * 0.5f;
	}

	// Enforce incompressibility and surface tension, double density relaxation
	// Algorithm 2 (4)
	for (size_t i = 0; i < n; i++) {
		float density = 0.0f, nearDensity = 0.0f;
		// Compute density and nearDensity
		for (size_t j = 0; j < n; j++) {	// for each neighboring pair ij
			if (j == i) continue;
			Vec2 rij = particles[j].pos - particles[i].pos;
			float dist2 = rij.squaredNorm();
			if (dist2 < radius * radius) {
				float dist = std::sqrt(dist2);
				float q = dist / radius;
				float omq = 1.0f - q;
				density += omq * omq;
				nearDensity += omq * omq * omq;
			}
		}
		// Compute pressure and nearPressure
		float pressure = stiffness * (density - restDensity);
		float nearPressure = stiffnessNear * nearDensity;
		if (pressure < 0) pressure = 0;

		Vec2 displacement(0.0f, 0.0f);
		for (size_t j = 0; j < n; j++) {
			if (j == i) continue;
			Vec2 rij = particles[j].pos - particles[i].pos;
			float dist2 = rij.squaredNorm();
			if (dist2 < radius * radius) {
				// Apply displacements
				float dist = std::sqrt(dist2);
				float q = dist / radius;
				float omq = 1.0f - q;
				Vec2 dir = (dist > 1e-6f) ? rij * (1.0f / dist) : Vec2(1.0f, 0.0f);
				Vec2 D = dir * ((pressure * omq + nearPressure * omq * omq) * dt * dt);
				particles[j].pos += D * 0.5f;
				displacement -= D * 0.5f;
			}
		}
		particles[i].pos += displacement;
	}
	
	// Recompute velocity from relaxed positions
	// Algorithm 1 (3)
	for (size_t i = 0; i < n; i++) {
		// Use previous position to compute next velocity
		particles[i].vel = (particles[i].pos - prevPos[i]) * (1.0f / dt);
	}
}