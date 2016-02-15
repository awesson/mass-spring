#include "SpringForce.h"
#include <GLUT/glut.h>

SpringForce::SpringForce(Particle* i_p1, Particle* i_p2, double i_dist, double i_ks, double i_kd) :
  p1(i_p1), p2(i_p2), dist(i_dist), ks(i_ks), kd(i_kd) {}

void SpringForce::draw()
{
  glBegin( GL_LINES );
  glColor3f(0.6, 0.7, 0.8);
  glVertex2f( p1->Position[0], p1->Position[1] );
  glColor3f(0.6, 0.7, 0.8);
  glVertex2f( p2->Position[0], p2->Position[1] );
  glEnd();
}

void SpringForce::add_force(){
        Vec2f dx = p1->Position - p2->Position;
        double norm_dx = norm(dx);
        double v_dx;
        // so that if the particles are on top of eachother the program does not blow up
        if(norm_dx == 0)
            v_dx = INF;
        else
            v_dx = (p1->Velocity - p2->Velocity) * dx / norm_dx;
        Vec2f force1 = -(ks*(norm_dx - dist) + kd*v_dx)*dx/norm_dx;
	p1->forces += force1;
        p2->forces -= force1;
}
