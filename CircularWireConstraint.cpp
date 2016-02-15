#include "CircularWireConstraint.h"
#include <GLUT/glut.h>

#define PI 3.1415926535897932384626433832795

static void draw_circle(const Vec2f & vect, float radius)
{
	glBegin(GL_LINE_LOOP);
	glColor3f(0.0,1.0,0.0); 
	for (int i=0; i<360; i=i+18)
	{
		float degInRad = i*PI/180;
		glVertex2f(vect[0]+cos(degInRad)*radius,vect[1]+sin(degInRad)*radius);
	}
	glEnd();
}

CircularWireConstraint::CircularWireConstraint(Particle *i_p, const Vec2f & i_center, const double i_radius) :
        p(i_p), center(i_center), radius(i_radius) {}

void CircularWireConstraint::draw()
{
        draw_circle(center, radius);
}

double CircularWireConstraint::get_C(){
    return norm2(p->Position - center) - pow(radius,2);
}

double CircularWireConstraint::get_Cdot(){
    return p->Velocity * (p->Position - center) / norm(p->Position - center);
}

Vec2f CircularWireConstraint::get_J(){
    Vec2f X = p->Position - center;
    if(norm(X) == 0)
        return Vec2f(INF, INF);
    return X / norm(X);
}

Vec2f CircularWireConstraint::get_Jdot(){
    Vec2f X = p->Position - center;
    if(norm(X) == 0)
        return Vec2f(INF, INF);
    return (p->Velocity - (p->Velocity * X) * X / norm2(X)) / norm(X);
}

double CircularWireConstraint::get_mass(){
    return p->mass;
}

int CircularWireConstraint::get_id(){
    return p->id;
}
