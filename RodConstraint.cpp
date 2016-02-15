#include "RodConstraint.h"
#include <GLUT/glut.h>

RodConstraint::RodConstraint(Particle* i_p1, Particle* i_p2, double i_dist) :
  p1(i_p1), p2(i_p2), dist(i_dist) {}

void RodConstraint::draw()
{
  glBegin( GL_LINES );
  glColor3f(0.8, 0.7, 0.6);
  glVertex2f( p1->Position[0], p1->Position[1] );
  glColor3f(0.8, 0.7, 0.6);
  glVertex2f( p2->Position[0], p2->Position[1] );
  glEnd();

}

double RodConstraint::get_C(){
    return norm(p1->Position - p2->Position) - dist;
}

double RodConstraint::get_Cdot(){
    Vec2f X = p1->Position - p2->Position;
    // so that if the particles are on top of eachother the program does not die
    if(norm(X) == 0)
        return INF;
    return (p1->Velocity - p2->Velocity) * X / norm(X);
}

// returns dC/dx1 since dC/dx2 = -dC/dx1
Vec2f RodConstraint::get_J(){
    Vec2f X = p1->Position - p2->Position;
    // so that if the particles are on top of eachother the program does not die
    if(norm(X) == 0)
        return Vec2f(INF, INF);
    return X / norm(X);
}

// returns d(dC/dx1)/dt since d(dC/dx2)/dt = -d(dC/dx1)/dt
Vec2f RodConstraint::get_Jdot(){
    Vec2f X = p1->Position - p2->Position;
    double norm_X = norm(X);
    Vec2f V = p1->Velocity - p2->Velocity;
    // so that if the particles are on top of eachother the program does not die
    if(norm_X == 0)
        return Vec2f(INF, INF);
    return V / norm_X - X * (X * V) / (norm_X * norm2(X));
}

double RodConstraint::get_mass1(){
    return p1->mass;
}

int RodConstraint::get_id1(){
    return p1->id;
}

double RodConstraint::get_mass2(){
    return p2->mass;
}

int RodConstraint::get_id2(){
    return p2->id;
}
