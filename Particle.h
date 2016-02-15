#pragma once

#include <gfx/vec2.h>

class Particle
{
public:

        Particle(const Vec2f & ConstructPos, double mass, int i_id);
	virtual ~Particle(void);

	void reset();
	void draw();

	Vec2f ConstructPos;
	Vec2f Position;
	Vec2f Velocity;
	Vec2f forces;
        // the derivitive of position and velocity respectivly
        Vec2f deriv_position;
        Vec2f deriv_velocity;
        // particle's numbering used in calculating lambda
        int id;
        double mass;
};
