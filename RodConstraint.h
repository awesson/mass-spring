#pragma once

#include "Particle.h"
#define INF 1e5

class RodConstraint {
 public:
  RodConstraint(Particle *i_p1, Particle* i_p2, double i_dist);

  void draw();
  double get_C();
  double get_Cdot();
  Vec2f get_J();
  Vec2f get_Jdot();
  double get_mass1();
  int get_id1();
  double get_mass2();
  int get_id2();

 private:

  Particle * const p1;
  Particle * const p2;
  double const dist;
};
