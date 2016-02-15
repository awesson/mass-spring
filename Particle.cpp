#include "Particle.h"
#include <GLUT/glut.h>

Particle::Particle(const Vec2f & i_ConstructPos, double i_mass, int i_id) :
	ConstructPos(i_ConstructPos), Position(Vec2f(0.0, 0.0)),
        Velocity(Vec2f(0.0, 0.0)), forces(Vec2f(0.0, 0.0)),
        deriv_position(Vec2f(0.0, 0.0)),
        deriv_velocity(Vec2f(0.0, 0.0)), id(i_id), mass(i_mass)
{
}

Particle::~Particle(void)
{
}

void Particle::reset()
{
	Position = ConstructPos;
	Velocity = Vec2f(0.0, 0.0);
        forces = Vec2f(0.0, 0.0);
}
void Particle::draw()
{
	const double h = 0.03;
	glColor3f(1.f, 1.f, 1.f); 
	glBegin(GL_QUADS);
	glVertex2f(Position[0]-h/2.0, Position[1]-h/2.0);
	glVertex2f(Position[0]+h/2.0, Position[1]-h/2.0);
	glVertex2f(Position[0]+h/2.0, Position[1]+h/2.0);
	glVertex2f(Position[0]-h/2.0, Position[1]+h/2.0);
	glEnd();
}
