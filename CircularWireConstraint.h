#pragma once

#include "Particle.h"
#define INF 1e5

class CircularWireConstraint {
 public:
  CircularWireConstraint(Particle *i_p, const Vec2f & i_center, const double i_radius);

  void draw();
  double get_C();
  double get_Cdot();
  Vec2f get_J();
  Vec2f get_Jdot();
  double get_mass();
  int get_id();

 private:

  Particle * const p;
  Vec2f const center;
  double const radius;
};
